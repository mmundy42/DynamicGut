from os import listdir
from os.path import splitext, join
import six.moves.cPickle as pickle
from warnings import warn
import cobra.io as io
from micom.util import load_model

_write_funcs = {
    '.xml': io.read_sbml_model,
    '.gz': io.read_sbml_model,
    '.mat': io.load_matlab_model,
    '.json': io.load_json_model,
    '.pickle': lambda fn: pickle.load(open(fn, 'rb'))
}


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
        _, ext = splitext(filename)
        if ext in _write_funcs:
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
        _, ext = splitext(name)
        _write_funcs[ext](model, name)
    return


def update_agora_models(model_file_names):

    # How to find biomass exchange reaction?
    for name in model_file_names:
        model = load_model(name)
        try:
            biomass_rxn = model.reactions.get_by_id('EX_biomass_LPAREN_e_RPAREN_')
            biomass_rxn.id = 'DM_biomass_LPAREN_e_RPAREN_'
            biomass_rxn.name = 'Demand Biomass c0'
            # cobra.io.write_sbml_model(model, name)
        except KeyError:
            warn('Model {0} in file {1} does not have biomass reaction'.format(model.id, name))
    return
