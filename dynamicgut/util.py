from os import listdir
from os.path import join
import json

from mminte import load_model_from_file, save_model_to_file


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
