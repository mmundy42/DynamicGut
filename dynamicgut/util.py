from os.path import join
import json
import pandas as pd
import re
from six import iteritems, iterkeys
from collections import defaultdict

from cobra.exceptions import OptimizationError
from cobra.io import write_sbml_model
from micom.util import load_model
from micom import Community, load_pickle

from .constants import pair_rate_columns, NO_GROWTH, ALMOST_ZERO
from .logger import logger


def check_for_growth(model_file_name, solver=None):
    """ Optimize a model and check for growth under conditions set in model.

    Parameters
    ----------
    model_file_name : str
        Path to input model file
    solver : str, optional
        Name of solver to use for optimization

    Returns
    -------
    dict
        Dictionary with summary of optimization results
    """

    model = load_model(model_file_name)
    if solver is not None:
        model.solver = solver
    summary = {'grows': True, 'message': None}
    try:
        value = model.slim_optimize()  # @todo use a timeout
        if value <= NO_GROWTH:
            summary['grows'] = False
            summary['message'] = 'Model {0} in file {1} does not produce growth under given conditions' \
                .format(model.id, model_file_name)
    except OptimizationError:
        summary['grows'] = False
        summary['message'] = 'Model {0} in file {1} fails to optimize'.format(model.id, model_file_name)
    return summary


def optimize_single_model(model_file_name, diet, compartment='e', solver=None):
    """ Optimize a single species model on a given diet conditions.

    This function can be used as a target function in a multiprocessing pool.

    Since the model is read from a file each time the function runs there is
    no need to revert the model after the optimization. This lowers the memory
    requirements of the simulation but there are more accesses of the file
    system.

    Parameters
    ----------
    model_file_name : str
        Path to single species model file
    diet : dict
        Dictionary with base metabolite ID as key and bound as value
    compartment : str
        Compartment ID for extracellular compartment in exchange reactions
    solver : str, optional
        Name of solver to use for optimization or None to use default solver

    Returns
    -------
    dict
        Dictionary with details on solution where exchange reaction fluxes are
        returned by base metabolite ID
    """

    # Load, apply the medium, and optimize the model.
    logger.debug('Optimizing single model "%s"', model_file_name)
    model = load_model(model_file_name)
    exchange_reactions = model.reactions.query(lambda x: x.startswith('EX_'), 'id')
    apply_medium(model, make_medium(exchange_reactions, diet, compartment))
    if solver is not None:
        model.solver = solver
    solution = model.optimize()  # @todo add timeout

    # Get the details on the solution.
    # @todo Need to check for timeout
    details = dict()
    details['model_id'] = model.id
    details['status'] = solution.status
    details['objective_value'] = solution.objective_value
    details['exchange_fluxes'] = dict()
    if solution.status == 'optimal':
        suffix = re.compile(r'_([{}])$'.format(compartment))
        for rxn in exchange_reactions:
            if solution.fluxes[rxn.id] != 0.0:
                met = next(iter(iterkeys(rxn.metabolites)))
                met_id = re.sub(suffix, '', met.id)
                details['exchange_fluxes'][met_id] = solution.fluxes[rxn.id]
    logger.debug('Model %s: %s %f', details['model_id'], details['status'], details['objective_value'])
    return details


def create_pair_model(pair, output_folder, to_sbml=False):
    """ Create a two species community model.

    Parameters
    ----------
    pair : tuple
        Each element is a path to a single species model file
    output_folder : str
        Path to output folder where community model JSON file is saved
    to_sbml : boolean, optional
        When True also save community model in SBML format

    Returns
    -------
    str
        Path to two species community model file
    """

    # Load the single species models just to get the model IDs.
    model_a = load_model(pair[0])
    model_b = load_model(pair[1])
    community_id = '{0}+{1}'.format(model_a.id, model_b.id)

    # Create a two species community model and save to a file.
    taxonomy = pd.DataFrame({"id": [model_a.id, model_b.id], "file": pair})
    community = Community(taxonomy, id=community_id, progress=False)
    community_file_name = join(output_folder, community.id + '.pickle')
    if to_sbml:
        write_sbml_model(community, join(output_folder, community.id + '.sbml'))
    community.to_pickle(community_file_name)
    return community_file_name


def map_single_to_pairs(pair_file_names):
    """ Find all of the pair community models that each single species model is a member of.

    Parameters
    ----------
    pair_file_names : list of str
        Each element is a path to a two species community model file

    Returns
    -------
    dict
        Dictionary keyed by single species model ID where value is list of paths
        to two species community model files that single species is a member of
    """

    single_to_pairs = defaultdict(list)
    for model_file_name in pair_file_names:
        pair_model = load_pickle(model_file_name)
        single_to_pairs[pair_model.taxonomy['id'][0]].append(model_file_name)
        single_to_pairs[pair_model.taxonomy['id'][1]].append(model_file_name)
    return single_to_pairs


def optimize_pair_model(model_file_name, diet, solver=None):
    """ Optimize a two species community model.

    This function can be used as a target function in a multiprocessing pool.

    Since the model is read from a file each time the function runs there is
    no need to revert the model after the optimization. This lowers the memory
    requirements of the simulation but there are more accesses of the file
    system.

    Current approach is to calculate the effect of species B on the growth of
    species A using the equation "G_ta / G_a" where G_ta is the growth rate of
    species A in the presence of species B and G_a is the growth rate of species
    A in the absence of species B. The same approach is used to calculate the
    effect of species A on the growth of species B.

    Note that an infeasible solution is considered the same as no growth.

    Parameters
    ----------
    model_file_name : str
        Path to two species community model file
    diet : dict
        Dictionary with base metabolite ID as key and bound as value
    solver : str, optional
        Name of solver to use for optimization or None to use default solver

    Returns
    -------
    pandas.Series
        Growth rate details for interaction between two species in pair
    """

    # Optimize the model the two species community model.
    logger.debug('Optimizing pair model "%s"', model_file_name)
    pair_model = load_pickle(model_file_name)
    a_id = pair_model.taxonomy['id'][0]
    b_id = pair_model.taxonomy['id'][1]
    pair_model.medium = make_medium(pair_model.exchanges, diet, 'm')
    if solver is not None:
        pair_model.solver = solver
    t_solution = pair_model.optimize(slim=True)
    a_alone = pair_model.optimize_single(a_id)  # @todo Check with Christian on slim
    b_alone = pair_model.optimize_single(b_id)

    # Round very small growth rates to zero.
    if a_alone < NO_GROWTH:
        a_alone = 0.
    if b_alone < NO_GROWTH:
        b_alone = 0.

    # Evaluate the interaction between the two species.
    if t_solution.status == 'optimal':
        a_together = t_solution.members.growth_rate[a_id]
        if a_together < NO_GROWTH:
            a_together = 0.
        b_together = t_solution.members.growth_rate[b_id]
        if b_together < NO_GROWTH:
            b_together = 0.
    else:
        a_together = 0.0
        b_together = 0.0

    if a_alone != 0.0:
        alone = a_alone
    else:
        alone = ALMOST_ZERO
    if a_together != 0.0 or a_alone != 0.0:
        a_effect = a_together / alone  # See note above for description
    else:
        a_effect = 0.0

    if b_alone != 0.0:
        alone = b_alone
    else:
        alone = ALMOST_ZERO
    if b_together != 0.0 or b_alone != 0.0:
        b_effect = b_together / alone  # See note above for description
    else:
        b_effect = 0.0

    logger.debug('Model %s: %f %f %f %f', pair_model.id, a_together, a_alone, b_together, b_alone)
    return pd.Series([a_id, b_id, a_together, a_alone, a_effect, b_together, b_alone, b_effect],
                     index=pair_rate_columns)


def make_medium(exchange_reactions, diet, compartment):
    """ Create a medium for a model from a diet.

    A global medium is a dictionary of bounds for metabolites which is converted
    to exchange reactions that can be applied to a model.

    Parameters
    ----------
    exchange_reactions : cobra.core.DictList
        List of exchange reactions from a model
    diet : dict
        Dictionary with base metabolite ID as key and bound as value
    compartment : str
        ID of compartment used as a suffix on exchange reaction IDs

    Returns
    -------
    dict
        Dictionary of the bounds for exchange reactions in pair community model
    """

    suffix = re.compile(r'_([{}])$'.format(compartment))
    medium = dict()
    for rxn in exchange_reactions:
        met = next(iter(iterkeys(rxn.metabolites)))
        met_id = re.sub(suffix, '', met.id)
        try:
            medium[rxn.id] = diet[met_id]
        except KeyError:
            pass
    return medium


def apply_medium(model, medium):
    """ Apply a medium to a model to set the metabolites that can be consumed.

    This function is adapted from the cobra.core.Model.medium setter in cobra 0.6
    with two differences: (1) if a reaction is in the medium but not in the
    model, the reaction is ignored (2) when turning off reactions in the model
    and not in the medium, only exchange reactions with the prefix `EX_` are
    considered (instead of all boundary reactions).

    Parameters
    ----------
    model : cobra.core.Model
        Model to apply medium to
    medium : dict
        Dictionary with exchange reaction ID as key and bound as value
    """

    def set_active_bound(reaction, bound):
        if reaction.reactants:
            reaction.lower_bound = -bound
        elif reaction.products:
            reaction.upper_bound = bound

    # Set the given media bounds.
    medium_reactions = set()
    for reaction_id, bound in iteritems(medium):
        try:
            reaction = model.reactions.get_by_id(reaction_id)
            medium_reactions.add(reaction)
            set_active_bound(reaction, bound)
        except KeyError:
            pass

    # The boundary attribute of a cobra.core.Reaction also includes demand and
    # sink reactions that we don't want turned off.
    exchange_reactions = set(model.reactions.query(lambda x: x.startswith('EX_'), 'id'))

    # Turn off reactions not present in medium.
    for reaction in (exchange_reactions - medium_reactions):
        set_active_bound(reaction, 0)

    return
