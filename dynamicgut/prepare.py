import six
from os import listdir, makedirs
from os.path import join, exists
from warnings import warn
from multiprocessing import Pool, cpu_count

from cobra.core import DictList
from cobra.io import save_json_model
from mminte import create_interaction_models, load_model_from_file


def prepare_simulation(model_file_names, single_model_folder, pair_model_folder, optimize=False, n_processes=None):
    """ Prepare single species models for simulation.

    Parameters
    ----------
    model_file_names : list of str
        Paths to input single species model files
    single_model_folder : str
        Path to folder for storing updated single species model files
    pair_model_folder : str
        Path to folder for storing two species community models
    optimize : bool
        When True, confirm single species models optimize successfully before creating community models
    n_processes: int, optional
        Number of processes in job pool

    Returns
    -------
    list of str
        List of paths to two species community models
    """

    # Create a job pool.
    if n_processes is None:
        n_processes = min(cpu_count(), 4)
    pool = Pool(n_processes)

    # Create folders for output models if needed.
    if not exists(single_model_folder):
        makedirs(single_model_folder)
    if not exists(pair_model_folder):
        makedirs(pair_model_folder)

    # Confirm all of the single species input models produce growth as provided.
    if optimize:
        result_list = [pool.apply_async(optimize_model, (file_name,))
                       for file_name in model_file_names]
        for result in result_list:
            summary = result.get()
            if summary['status'] != 'optimal':
                warn('Model {0} in file {1} fails to optimize'.format(summary['id'], summary['file_name']))
            elif summary['f'] <= 0.0001:
                warn('Model {0} in file {1} does not produce growth under given conditions'
                     .format(summary['id'], summary['file_name']))

    # Extract all of the exchange reactions from the single species models.
    # @todo Could also save the exchange reactions in a separate model.
    exchange_reactions = get_exchange_reactions(model_file_names)

    # Make sure each single species model has a complete set of exchange reactions. Then each
    # two species community model will have a complete set of exchange reactions. Updated
    # single species models are stored in the specified folder.
    result_list = [pool.apply_async(add_exchange_reactions, (file_name, exchange_reactions, single_model_folder))
                   for file_name in model_file_names]
    single_file_names = [result.get() for result in result_list]
    pool.close()

    # Create all of the pair models and store in specified folder.
    create_interaction_models(single_file_names, pair_model_folder, n_processes=n_processes)
    return


def find_models_in_folder(source_folder):
    """ Get a list of model file path names based on supported extensions.

        Note that the contents of a file might not be a model so this should
        only be used for folders where all of the files are model files.

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


def optimize_model(model_file_name):
    """ Optimize a model and return summary results.

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
    return {
        'file_name': model_file_name,
        'id': model.id,
        'status': solution.status,
        'f': solution.f
    }


def add_exchange_reactions(model_file_name, exchange_reactions, model_folder):
    """ Add missing exchange reactions to a model.

    Parameters
    ----------
    model_file_name : str
        Path to input model file
    exchange_reactions : cobra.core.DictList
        List of all exchange reactions
    model_folder : str
        Path to folder for storing updated model file

    Returns
    -------
    str
        Path to updated model file
    """

    # Load the input model file.
    model = load_model_from_file(model_file_name)

    # Set the bounds of the exchange reaction to the default value or add
    # a missing exchange reaction to the model.
    for reaction in exchange_reactions:
        try:
            model.reactions.get_by_id(reaction.id).bounds = (-1000., 1000.)
        except KeyError:
            model.add_reaction(reaction)

    # Updated model is stored in JSON format.
    json_file_name = join(model_folder, model.id+'.json')
    save_json_model(model, json_file_name)

    return json_file_name


def get_exchange_reactions(model_file_names):
    """ Get a list of all unique exchange reactions in a set of models.

    Parameters
    ----------
    model_file_names : list of str
        Paths to model files

    Returns
    -------
    cobra.core.DictList
        List of cobra.Reaction objects for exchange reactions in input models
    """

    all_exchange_reactions = DictList()

    for file_name in model_file_names:
        model = load_model_from_file(file_name)

        # Get the exchange reactions from the model.
        exchange_reactions = model.reactions.query(lambda x: x.startswith('EX_'))
        for reaction in exchange_reactions:
            if not all_exchange_reactions.has_id(reaction.id):
                all_exchange_reactions.append(copy_exchange_reaction(reaction))

    return all_exchange_reactions


def copy_exchange_reaction(reaction):
    """ Make a copy of an exchange reaction from an exchange reaction in a source model.

    Parameters
    ----------
    reaction : cobra.Reaction
        Source exchange reaction

    Returns
    -------
    cobra.Reaction
        Copy of exchange reaction with default values
    """

    # Confirm reaction has only one metabolite and metabolite coefficient is negative.
    if len(reaction.metabolites) > 1:
        warn('Model {0} exchange reaction {1} has {2} metabolites (expected one metabolite)'
             .format(reaction.model.id, reaction.id, len(reaction.metabolites)))
    for metabolite in reaction.metabolites:
        if reaction.metabolites[metabolite] >= 0:
            warn('Model {0} exchange reaction {1} metabolite {2} has positive coefficient (expected negative)'
                 .format(reaction.model.id, reaction.id, metabolite.id))

    # Make a copy of the reaction. Since exchange reactions are shared by all species in the
    # community model, set the lower and upper bounds to default values.
    copy_reaction = reaction.copy()
    copy_reaction.lower_bound = -1000.
    copy_reaction.upper_bound = 1000.
    copy_reaction.objective_coefficient = 0.

    # Confirm metabolite is in extracellular compartment.
    metabolite = six.next(six.iterkeys(copy_reaction.metabolites))
    if metabolite.compartment != 'e':
        warn('Model {0} exchange reaction {1} metabolite {2} has wrong compartment {3} (expected "e")'
             .format(reaction.model.id, reaction.id, metabolite.id, metabolite.compartment))

    return copy_reaction
