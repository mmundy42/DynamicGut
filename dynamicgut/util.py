from os import listdir
from os.path import join

from mminte import load_model_from_file


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
    """ Get the set of unique exchange reaction IDs from a list of models.

    Parameters
    ----------
    model_file_names : list of str
        List of path names to model files

    Returns
    -------
    set
        Set of exchange reaction IDs in input models
    """

    all_exchange_reactions = set()
    model_list = [load_model_from_file(name) for name in model_file_names]
    for model in model_list:
        exchange_reactions = model.reactions.query(lambda x: x.startswith('EX_'), 'id')
        for reaction in exchange_reactions:
            all_exchange_reactions.add(reaction.id)
    return all_exchange_reactions
