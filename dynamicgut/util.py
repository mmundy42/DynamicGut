from os import listdir
from os.path import join
import json
import logging

from cobra.util.solver import linear_reaction_coefficients, set_objective

from mminte import load_model_from_file, save_model_to_file, apply_medium

# Logger for this module
LOGGER = logging.getLogger(__name__)


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

    source_models = list()
    for filename in listdir(source_folder):
        if filename.endswith('.mat') or filename.endswith('.xml') or \
                filename.endswith('.sbml') or filename.endswith('.json'):
            source_models.append(join(source_folder, filename))
    return source_models


def check_for_growth(model_file_name):
    """ Optimize a model and check for growth under conditions set in model.

    Parameters
    ----------
    model_file_name : str
        Path to input model file

    Returns
    -------
    dict
        Dictionary with summary of optimization results
    """

    model = load_model_from_file(model_file_name)
    solution = model.optimize()
    summary = {'grows': True, 'message': None}
    if solution.status != 'optimal':
        summary['grows'] = False
        summary['message'] = 'Model {0} in file {1} fails to optimize'.format(model.id, model_file_name)
    elif solution.objective_value <= 0.0001:
        summary['grows'] = False
        summary['message'] = 'Model {0} in file {1} does not produce growth under given conditions' \
                             .format(model.id, model_file_name)
    return summary


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
        model = load_model_from_file(name)
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
        model = load_model_from_file(name)
        exchange_reactions = model.reactions.query(lambda x: x.startswith('EX_'), 'id')
        for rxn in exchange_reactions:
            if rxn.id not in diet:
                if bound is None:
                    diet[rxn.id] = get_active_bound(rxn)
                else:
                    diet[rxn.id] = bound
    json.dump(diet, open(diet_file_name, 'w'), indent=4)
    return


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
        model = load_model_from_file(name)
        model.id = prefix + model.id
        save_model_to_file(model, name)
    return


def optimize_for_species(model_file_name, species_id, medium, time_point_folder):
    """ Knock out the other species from a two species model, optimize, and save results.
    """

    LOGGER.info('Loading model {0} to optimize for {1}'.format(model_file_name, species_id))
    pair_model = load_model_from_file(model_file_name)

    # Figure out the species to knock out in the pair community model.
    species_index = -1
    for index in range(len(pair_model.notes['species'])):
        if pair_model.notes['species'][index]['id'] == species_id:
            species_index = index
    if species_index < 0:
        raise Exception('Species {0} is not a member of the community'.format(species_id))
    if species_index == 0:
        knockout_index = 1
    else:
        knockout_index = 0
    knockout_id = pair_model.notes['species'][knockout_index]['id']
    LOGGER.info('Going to knock out {0} from index {1}'.format(knockout_id, knockout_index))

    with pair_model:
        # Apply the medium.
        apply_medium(pair_model, medium)

        # Knock out all of the reactions for the specified species.
        knockout_reactions = pair_model.reactions.query(lambda r: r.startswith(knockout_id), 'id')
        for reaction in knockout_reactions:
            reaction.knock_out()

        # Remove the species objective from the community model objective.
        knockout_objective = pair_model.reactions.get_by_id(pair_model.notes['species'][knockout_index]['objective'])
        linear_coefficients = linear_reaction_coefficients(pair_model)
        del linear_coefficients[knockout_objective]
        set_objective(pair_model, linear_coefficients)
        save_model_to_file(pair_model, join(time_point_folder, pair_model.id+'.json'))

        # Optimize the community model with the specified species knocked out.
        solution = pair_model.optimize()
        solution.fluxes.to_json(join(time_point_folder, pair_model.id+'-solution.json'))  # , orient='records', lines=True

    return solution
