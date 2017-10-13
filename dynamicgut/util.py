from os import listdir
from os.path import join
import json
import logging
import pandas as pd
import re
from six import iteritems
from collections import defaultdict

from cobra.exceptions import OptimizationError
from micom.util import load_model
from micom import Community, load_pickle

from .constants import pair_rate_columns, NO_GROWTH, ALMOST_ZERO

# Logger for this module
LOGGER = logging.getLogger(__name__)


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
        value = model.slim_optimize()
        if value <= NO_GROWTH:
            summary['grows'] = False
            summary['message'] = 'Model {0} in file {1} does not produce growth under given conditions' \
                .format(model.id, model_file_name)
    except OptimizationError:
        summary['grows'] = False
        summary['message'] = 'Model {0} in file {1} fails to optimize'.format(model.id, model_file_name)
    return summary


def optimize_single_model(model_file_name, medium, compartment='e', solver=None):
    """ Optimize a single species model on a given medium.

    This function can be used as a target function in a multiprocessing pool.

    Note that we chose to read the model from a file each time instead of loading
    the model into memory once at the beginning of the simulation. This lowers
    the memory requirements of the simulation and there is no need to revert the
    model after the optimization. But there are more accesses of the file system.

    Parameters
    ----------
    model_file_name : str
        Path to single species model file
    medium : dict
        Dictionary with global metabolite ID as key and bound as value
    compartment : str
        Compartment ID for extracellular compartment in exchange reactions
    solver : str, optional
        Name of solver to use for optimization or None to use default solver

    Returns
    -------
    dict
        Dictionary with details on solution where exchange reaction fluxes are
        returned by global metabolite ID
    """

    # Load, apply the medium, and optimize the model.
    model = load_model(model_file_name)
    exchange_reactions = model.reactions.query(lambda x: x.startswith('EX_'), 'id')
    apply_medium(model, make_medium(exchange_reactions, medium, compartment))
    if solver is not None:
        model.solver = solver
    solution = model.optimize()

    # Get the details on the solution.
    details = dict()
    details['model_id'] = model.id
    details['status'] = solution.status
    details['objective_value'] = solution.objective_value
    details['exchange_fluxes'] = dict()
    if solution.status == 'optimal':
        suffix = re.compile(r'_([{}])$'.format(compartment))
        for rxn in exchange_reactions:
            if solution.fluxes[rxn.id] != 0.0:
                met_id = re.sub(suffix, '', rxn.metabolites.keys()[0].id)
                details['exchange_fluxes'][met_id] = solution.fluxes[rxn.id]
    return details


def create_pair_model(pair, output_folder):
    """ Create a two species community model.

    Parameters
    ----------
    pair : tuple
        Each element is a path to a single species model file
    output_folder : str
        Path to output folder where community model JSON file is saved

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
    community = Community(taxonomy, id=community_id)
    community_file_name = join(output_folder, community.id + '.pickle')
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


def optimize_pair_model(model_file_name, medium, solver=None):
    """ Optimize a two species community model.

    This function can be used as a target function in a multiprocessing pool.
    Since the model is read from a file each time the function runs there is
    no need to revert the model after the optimization.

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
    medium : dict
        Dictionary with global metabolite ID as key and bound as value
    solver : str, optional
        Name of solver to use for optimization or None to use default solver

    Returns
    -------
    pandas.Series
        Growth rate details for interaction between two species in pair
    """

    # Optimize the model the two species community model.
    pair_model = load_pickle(model_file_name)
    a_id = pair_model.taxonomy['id'][0]
    b_id = pair_model.taxonomy['id'][1]
    pair_model.medium = make_medium(pair_model.exchanges, medium, 'm')
    if solver is not None:
        pair_model.solver = solver
    t_solution = pair_model.optimize(slim=True)
    a_alone = pair_model.optimize_single(a_id)
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

    return pd.Series([a_id, b_id, a_together, a_alone, a_effect, b_together, b_alone, b_effect],
                     index=pair_rate_columns)


def get_exchange_reaction_ids(model_file_names):
    """ Get the set of unique exchange reaction IDs and model IDs from a list of models.

    Parameters
    ----------
    model_file_names : list of str
        List of path names to model files

    Returns
    -------
    set
        Set of exchange reaction IDs in input models
    list
        List of model IDs from input models
    """

    all_exchange_reactions = set()
    model_ids = list()
    for name in model_file_names:
        model = load_model(name)
        model_ids.append(model.id)
        exchange_reactions = model.reactions.query(lambda x: x.startswith('EX_'), 'id')
        for rxn in exchange_reactions:
            all_exchange_reactions.add(rxn.id)
    return all_exchange_reactions, model_ids


def make_diet_from_models(model_file_names, diet_file_name, bound=None):
    """ Make a diet file from the exchange reactions in a list of models.

    Parameters
    ----------
    model_file_names : list of str
        List of path names to model files
    diet_file_name : str
        Path to file to store diet conditions in JSON format
    bound : float, optional
        Bound to set on every exchange reaction, when None use bound from first
        model that contains exchange reaction
    """

    def get_active_bound(reaction):
        """ For an active boundary reaction, return the relevant bound. """
        if reaction.reactants:
            return -reaction.lower_bound
        elif reaction.products:
            return reaction.upper_bound

    if bound is not None:
        bound = float(abs(bound))
    diet = dict()
    for name in model_file_names:
        model = load_model(name)
        exchange_reactions = model.reactions.query(lambda x: x.startswith('EX_'), 'id')
        for rxn in exchange_reactions:
            if rxn.id not in diet:
                if bound is None:
                    diet[rxn.id] = get_active_bound(rxn)
                else:
                    diet[rxn.id] = bound
    json.dump(diet, open(diet_file_name, 'w'), indent=4)
    return


def make_medium(exchange_reactions, global_medium, compartment):
    """ Create a medium for a model from a global medium.

    A global medium is a dictionary of bounds for metabolites which is converted
    to exchange reactions that can be applied to a model.

    Parameters
    ----------
    exchange_reactions : cobra.core.DictList
        List of exchange reactions from a model
    global_medium : dict
        Dictionary of the bounds for each metabolite in an exchange reaction
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
        met_id = re.sub(suffix, '', rxn.metabolites.keys()[0].id)
        try:
            medium[rxn.id] = global_medium[met_id]
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


def find_models_in_folder(source_folder):
    """ Get a list of model file path names based on supported extensions.

    Note that the contents of a file with one of the supported extensions might
    not be a model so this should only be used for folders where all of the files
    are model files.

    Parameters
    ----------
    source_folder : str
        Path to folder with source model files

    Returns
    -------
    list of str
        List of path names to source model files
    """

    # @todo Could use write_funcs dict here
    source_models = list()
    for filename in listdir(source_folder):
        if filename.endswith('.mat') or filename.endswith('.xml') or \
                filename.endswith('.sbml') or filename.endswith('.json'):
            source_models.append(join(source_folder, filename))
    return source_models


def set_model_id_prefix(model_file_names, prefix='M'):
    """ Set a prefix on model IDs for all models in a list of models.

    Model IDs must start with an alphabetic character so they are interpreted as
    strings in data frames. Models created by ModelSEED typically use a PATRIC
    genome ID as the model ID which is a number.

    Parameters
    ----------
    model_file_names : list of str
        List of path names to model files
    prefix : str, optional
        String to use as prefix for model IDs
    """

    for name in model_file_names:
        model = load_model(name)
        model.id = prefix + model.id
        # @todo Need to handle saving models
        # save_model_to_file(model, name)
    return
