from os import listdir
from os.path import splitext, join
import six.moves.cPickle as pickle
from warnings import warn
import json
from six import iterkeys
import re

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


def get_exchange_metabolite_ids(model_file_names, compartment='e'):
    """ Get the set of unique metabolite IDs from exchange reactions and the
        model IDs from a list of models.

    Parameters
    ----------
    model_file_names : list of str
        List of path names to model files
    compartment : str
        Compartment ID for extracellular compartment in exchange reactions

    Returns
    -------
    set
        Set of base metabolite IDs used in exchange reactions of input models
    list
        List of model IDs from input models
    """

    suffix = re.compile(r'_([{}])$'.format(compartment))
    all_exchange_metabolites = set()
    model_ids = list()
    for name in model_file_names:
        model = load_model(name)
        model_ids.append(model.id)
        exchange_reactions = model.reactions.query(lambda x: x.startswith('EX_'), 'id')
        for rxn in exchange_reactions:
            met = next(iter(iterkeys(rxn.metabolites)))
            all_exchange_metabolites.add(re.sub(suffix, '', met.id))
    return all_exchange_metabolites, model_ids


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
