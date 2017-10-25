from os import makedirs
from os.path import join, exists
from warnings import warn
from multiprocessing import Pool
import pandas as pd
import numpy as np
import json
from collections import defaultdict
from itertools import combinations
from six import iterkeys

from .micomutil import check_for_growth, create_pair_model, optimize_single_model, optimize_pair_model
from .modelutil import get_exchange_metabolite_ids
from .constants import single_rate_columns, pair_rate_columns, density_columns, NO_GROWTH
from .logger import logger

# For debugging pandas
# set_option('display.width', 1000)


def prepare(single_file_names, pair_model_folder, optimize=False, n_processes=None, solver=None, to_sbml=False):
    """ Prepare for simulation by creating two species community models.

    Parameters
    ----------
    single_file_names : list of str
        Paths to single species model files of community members
    pair_model_folder : str
        Path to folder for storing two species community models
    optimize : bool, optional
        When True, confirm single species models optimize successfully before creating community models
    n_processes: int, optional
        Number of processes in job pool or None to run without multiprocessing
    solver : str, optional
        Name of solver to use for optimization or None to use default solver (most likely glpk)
    to_sbml : boolean, optional
        When True also save community model in SBML format

    Returns
    -------
    list of str
        List of paths to two species community models
    """

    # If requested, confirm all of the single species input models produce growth as provided.
    if optimize:
        logger.info('Started checking %d single models for growth', len(single_file_names))
        if n_processes is None:
            summary_list = [check_for_growth(file_name, solver)
                            for file_name in single_file_names]
        else:
            pool = Pool(n_processes)
            result_list = [pool.apply_async(check_for_growth, (file_name, solver))
                           for file_name in single_file_names]
            summary_list = [result.get() for result in result_list]
            pool.close()
        for summary in summary_list:
            if not summary['grows']:
                warn(summary['message'])
        logger.info('Finished checking single models for growth')

    # Create folder for output pair models if needed.
    if not exists(pair_model_folder):
        makedirs(pair_model_folder)

    # Create all of the pair models and store in specified folder.
    logger.info('Started creating two species community models')
    pair_list = [pair for pair in combinations(single_file_names, 2)]
    if n_processes is None:
        pair_file_names = [create_pair_model(pair, pair_model_folder, solver, to_sbml)
                           for pair in pair_list]
    else:
        pool = Pool(n_processes)
        result_list = [pool.apply_async(create_pair_model, (pair, pair_model_folder, solver, to_sbml))
                       for pair in pair_list]
        pair_file_names = [result.get() for result in result_list]
        pool.close()
    logger.info('Finished creating %d two species community models', len(pair_file_names))
    return pair_file_names


def run_simulation(time_interval, single_file_names, pair_file_names, diet_file_name,
                   density_file_name, data_folder, time_step=0.5, k=1,
                   n_processes=None, solver=None):
    """ Run a simulation over a time interval.

    Parameters
    ----------
    time_interval : range
        Range of time points for running the simulation
    single_file_names : list of str
        List of path names to single species models
    pair_file_names : list of str
        List of path names to two species community model files
    diet_file_name : str
        Path to file with initial diet conditions in JSON format
    density_file_name : str
        Path to file with initial population densities in CSV format
    data_folder : str
        Path to folder for storing data generated at each time point
    time_step : float, optional
        Adjustment to time point where 1 is one hour, 0.5 is 30 minutes, etc.
    k : int, optional
        Maximum size of the population that the environment has the capacity to support
    n_processes: int, optional
        Number of processes in job pool or None to run without multiprocessing
    solver : str, optional
        Name of solver to use for optimization or None to use default solver
    """

    # Basic validation of input parameters.
    if time_step <= 0.0 or time_step > 1.0:
        raise ValueError('The value of the time_step parameter must be greater than 0 and less than or equal to 1')
    if len(single_file_names) < 2:
        raise ValueError('There must be at least two models in list of single species model files')
    if len(pair_file_names) < 1:
        raise ValueError('There must be at least one model in list of two species community model files')

    # Create the data folder if needed (exception is raised if there is a problem).
    if not exists(data_folder):
        makedirs(data_folder)

    # Get the initial population density values.
    logger.info('Reading initial population density from "%s"', density_file_name)
    density = pd.read_csv(density_file_name, dtype={'MODEL_ID': str, 'DENSITY': float})
    if not set(density_columns).issubset(density.columns):
        raise ValueError('Required columns {0} not available in initial population density file "{1}"'
                         .format(density_columns, density_file_name))
    invalid_fields = density.isnull().values.sum()
    if invalid_fields > 0:
        raise ValueError('There are {0} fields with invalid values in initial population density file "{1}"'
                         .format(invalid_fields, density_file_name))

    # Set diet for first time step by adding exchange reactions from the single
    # species models that are not in the initial diet. This allows metabolites
    # produced by a species to become available in the diet conditions during the
    # simulation.
    logger.info('Reading initial diet from "%s"', diet_file_name)
    diet = json.load(open(diet_file_name))
    initial_exchanges = set(iterkeys(diet))
    logger.info('Getting metabolites from exchange reactions in %d single species models', len(single_file_names))
    model_exchanges, model_ids = get_exchange_metabolite_ids(single_file_names)
    if initial_exchanges > model_exchanges:
        warn('Diet file "{0}" contains more exchange reactions than there are in single species models'
             .format(diet_file_name))
    missing_metabolites = model_exchanges - initial_exchanges
    logger.info('Adding %d missing metabolites to initial diet', len(missing_metabolites))
    for met_id in missing_metabolites:
        diet[met_id] = 0.0
        logger.debug('Added metabolite %s to initial diet', met_id)
    json.dump(diet, open(join(data_folder, 'initial-diet.json'), 'w'), indent=4)

    # Confirm the model IDs in the initial density file match the model IDs in the
    # list of single species models.
    if density.shape[0] != len(model_ids):
        for index, row in density.iterrows():
            if row['MODEL_ID'] not in model_ids:
                logger.error('Model ID "{0}" on line {1} of initial population density file "{2}" is not available '
                             'in list of single species models'.format(row['ID'], index + 2, density_file_name))
        for model_id in model_ids:
            if model_id not in density.ID.values:
                logger.error('Model ID "{0}" from list of single species models is not available in '
                             'initial population density file'.format(model_id))
        raise ValueError('Number of species ({0}) in initial density file does not match '
                         'number of single species models ({1})'.format(density.shape[0], len(model_ids)))
    if density.loc[density['MODEL_ID'].isin(model_ids)].shape[0] != len(model_ids):
        for index, row in density.iterrows():
            if row['MODEL_ID'] not in model_ids:
                logger.error('Model ID "{0}" on line {1} of initial population density file "{2}" is not available '
                             'in list of single species models'.format(row['ID'], index + 2, density_file_name))
        raise ValueError('One or more model IDs in initial density file do not match '
                         'model IDs in single species models')

    # Run the simulation over the specified time interval.
    for time_point in time_interval:
        # Start this time point.
        time_point_id = '{0:04d}'.format(time_point + 1)
        logger.info('[%s] STARTED TIME POINT', time_point_id)
        time_point_folder = join(data_folder, 'timepoint-' + time_point_id)
        if not exists(time_point_folder):
            makedirs(time_point_folder)
        diet_file_name = join(time_point_folder, 'diet-{0}.json'.format(time_point_id))
        json.dump(diet, open(diet_file_name, 'w'), indent=4)
        pair_rate_file_name = join(time_point_folder, 'pair-rates-{0}.csv'.format(time_point_id))
        effects_matrix_file_name = join(time_point_folder, 'effects-matrix-{0}.csv'.format(time_point_id))
        density_file_name = join(time_point_folder, 'density-{0}.csv'.format(time_point_id))
        single_rate_file_name = join(time_point_folder, 'single-rates-{0}.csv').format(time_point_id)
        next_diet_file_name = join(time_point_folder, 'next-diet-{0}.json'.format(time_point_id))

        # Calculate the growth rates for each two species model under the current diet conditions.
        growth_rates, alone = calculate_growth_rates(pair_file_names, diet, pair_rate_file_name,
                                                     time_point_id, n_processes, solver)

        # Create the effects matrix.
        effects_matrix = create_effects_matrix(growth_rates, effects_matrix_file_name, time_point_id)

        # Run Leslie-Gower algorithm to calculate new population densities.
        density = leslie_gower(effects_matrix, density, density_file_name, time_point_id, alone, k, time_step)

        # Get the exchange reaction fluxes from optimizing single species models.
        exchange_fluxes = get_exchange_fluxes(single_file_names, diet, single_rate_file_name,
                                              time_point_id, n_processes, solver)

        # Create diet conditions for next time point.
        diet = create_next_diet(diet, exchange_fluxes, density, next_diet_file_name, time_step, time_point_id)

    # Cleanup and store results from last time step in data folder.
    json.dump(diet, open(join(data_folder, 'final-diet.json'), 'w'), indent=4)
    density.to_csv(join(data_folder, 'final-density.csv'))

    return


def calculate_growth_rates(pair_file_names, current_diet, pair_rate_file_name, time_point_id,
                           n_processes=None, solver=None):
    """ Optimize the two species community models to get growth rates.

    Parameters
    ----------
    pair_file_names : list of str
        Paths to two species community model files
    current_diet : dict
        Dictionary with global metabolite ID as key and bound as value
    pair_rate_file_name : str
        Path to file for storing growth rates of two species community models
    time_point_id : str
        ID of current time point
    n_processes : int, optional
        Number of processes in job pool or None to run without multiprocessing
    solver : str, optional
        Name of solver to use for optimization or None to use default solver

    Returns
    -------
    pandas.DataFrame
        Growth rates
    dict
        Single species growth rates
    """

    logger.info('[%s] Optimizing %d two species community models for interactions',
                time_point_id, len(pair_file_names))
    if solver is not None:
        logger.debug('[%s] Solver %s ignored with micom community models', time_point_id, solver)

    # Optimize all of the two species community models on the current diet conditions.
    if n_processes is None:
        rate_list = [optimize_pair_model(file_name, current_diet)
                     for file_name in pair_file_names]
    else:
        pool = Pool(n_processes)
        result_list = [pool.apply_async(optimize_pair_model, (file_name, current_diet))
                       for file_name in pair_file_names]
        rate_list = [result.get() for result in result_list]
        pool.close()

    # Build a DataFrame with the pair growth rates.
    pair_rate = pd.DataFrame(columns=pair_rate_columns)
    for rate in rate_list:
        if rate['STATUS'] != 'optimal':
            logger.warn('[%s] Pair model %s+%s did not grow on current diet conditions, status is %s',
                        time_point_id, rate['A_ID'], rate['B_ID'], rate['STATUS'])
        pair_rate = pair_rate.append(rate, ignore_index=True)
    pair_rate.to_csv(pair_rate_file_name, index=False)

    # Build a dictionary with single species growth rates and confirm that the
    # values are consistent.
    alone_rate = dict()
    inconsistent = set()
    for index, row in pair_rate.iterrows():
        if row['A_ID'] in alone_rate:
            if not np.isclose(row['A_ALONE'], alone_rate[row['A_ID']]):
                inconsistent.add(row['A_ID'])
                # logger.debug('[{0}] Model {1} has inconsistent growth rates: {2} vs {3}'
                #              .format(time_point_id, row['A_ID'], row['A_ALONE'], alone_rate[row['A_ID']]))
        else:
            alone_rate[row['A_ID']] = row['A_ALONE']
        if row['B_ID'] in alone_rate:
            if not np.isclose(row['B_ALONE'], alone_rate[row['B_ID']]):
                inconsistent.add(row['B_ID'])
                # logger.debug('[{0}] Model {1} has inconsistent growth rates: {2} vs {3}'
                #              .format(time_point_id, row['B_ID'], row['B_ALONE'], alone_rate[row['B_ID']]))
        else:
            alone_rate[row['B_ID']] = row['B_ALONE']

    return pair_rate, alone_rate


def create_effects_matrix(pair_rate, effects_matrix_file_name, time_point_id):
    """ Create an effects matrix from growth rate data frame of all pairs.

    Each cell in an effects matrix is the effect on the growth of one species in
    the presence of another species. A row gives the magnitude  in growth of
    a specific species in the presence of all of the other species in the community.
    The diagonal in the matrix is always 1.

    Parameters
    ----------
    pair_rate : pandas.DataFrame
        Growth rate details for all pairs in the community
    effects_matrix_file_name : str
        Path to file for storing effects matrix
    time_point_id : str
        ID of current time point

    Returns
    -------
    pandas.DataFrame
        Effect of one species on the growth of another species
    """

    logger.info('[%s] Creating effects matrix', time_point_id)

    # Extract the species model IDs and effect values from the input data frame.
    a_id = pair_rate['A_ID'].tolist()
    b_id = pair_rate['B_ID'].tolist()
    a_effect = pair_rate['A_EFFECT'].tolist()
    b_effect = pair_rate['B_EFFECT'].tolist()

    # Build the output data frame.
    raw = pd.DataFrame({'EFFECT': a_id, 'BECAUSE_OF': b_id, 'CHANGE': a_effect})
    raw = raw.append(pd.DataFrame({'EFFECT': b_id, 'BECAUSE_OF': a_id, 'CHANGE': b_effect}),
                     ignore_index=True)
    effects = raw.pivot_table(index='EFFECT', columns='BECAUSE_OF', values='CHANGE')
    effects = effects.replace(np.nan, 1.0, regex=True)
    effects = effects.replace(0.0, 1e-25, regex=True)
    effects.to_csv(effects_matrix_file_name)
    return effects


def leslie_gower(effects_matrix, current_density, density_file_name, time_point_id, alone_rate, k=1, time_step=0.5):
    """ Run Leslie-Gower algorithm to update population density of each species.

    Leslie Equation 3-4 defines the population size of a species in a pair to be:

    N1(t+1) = lambda1 * N1(t) / 1 + alpha1 * N1(t) + gamma1 * N2(t)

    Terms are defined in the code below.

    Parameters
    ----------
    effects_matrix : pandas.DataFrame
         Effect of one species on the growth of another species
    current_density : pandas.DataFrame
        Current population densities for species in community
    density_file_name : str
        Path to file for storing population densities at this time point
    time_point_id : str
        ID of current time point
    alone_rate : dict
        Single species growth rates
    k : int, optional
        Maximum size of the population that the environment has the capacity to support
    time_step : float
        Size of time interval where the value is greater than 0 and less than or
        equal to 1 and 1 means one hour, 0.5 means 30 minutes, etc.

    Returns
    -------
    pandas.DataFrame
        Updated population densities for species in community
    """

    logger.info('[%s] Calculating population densities', time_point_id)

    # Confirm that the order of the species in the effects matrix matches the
    # order in the density data frame.
    if not np.array_equal(effects_matrix.axes[0].values, current_density['MODEL_ID'].values):
        # @todo Figure out a way to include the mismatched data in the exception or log it
        raise ValueError('Order of species in effects matrix does not match density data frame')

    # The "N1(t)" term in equation above is the population size at the beginning of the
    # time point and comes from the DENSITY column of the input data frame. Create a
    # 1-D array from the data frame.
    density = current_density['DENSITY'].values

    # Convert the dictionary of single species growth rates at this time point to
    # a 1-D array in the same order as the current density array.
    growth = np.array([alone_rate[row['MODEL_ID']] for index, row in current_density.iterrows()], dtype=float)

    # Growth rate is reported in units of one hour. Adjust the growth rate by the
    # size of the time step (e.g. reduce by half for a 30 minute time step).
    growth = growth * time_step

    # The "lambda1" term in the equation above is defined as "1 + G(t)" where G(t)
    # is the growth rate at time point t. Leslie Equation 2-3 defines lambda as
    # "1 + B(t) - D(t)" where B(t) is the birth rate at time point t and D(t) is
    # the death rate at time point t. DynamicGut assumes "G(t) = B(t) - D(t)".
    lambdas = 1 + growth

    # The "alpha1" term in the equation above is defined as "(lambda - 1) / K" where
    # K is the size of the population that the environment has the capacity to
    # support (Equation 3-3). DynamicGut assumes that the carrying capacity is the
    # same for all species in the microbial community.
    alphas = (lambdas - 1) / k

    # The "gamma1" term in equation above is the magnitude of the effect which each
    # species has on the rate of increase of the other species. Create a 2-D matrix
    # from the input data frame.
    effects = effects_matrix.values

    # Actual Effects (21 January 2017). This is the decrease in growth of the focal species due to the others.
    # Previously, we had the total growth of the focal species in the presence of the other, which i don't
    # think made much sense for being in the numerators.
    # Actually, I'm still not completely sure about this.
    # @todo Can this be removed because of statement below
    # @todo Can this be done when calculating effects matrix?
    # @todo An intial effect greater 1 turns into a negative number which decreases value of denominator
    # effects = 1 - effects
    effects = -(np.log10(effects))  # @TODO: Added by Lena (18 July 2017)
    effects = effects * alphas  # @TODO: Added by Lena (18 July 2017)

    # Calculate the "gamma1 * N2(t)" term in the equation above. The resulting
    # total_effects vector represents the total effect of all other species on
    # the focal species. The matrix multiplication of the 2-D effects matrix and
    # the 1-D density array, creates a 1-D array.
    total_effects = np.dot(effects, density)

    # Calculate the population size at the end of the time point using the equation
    # given above.
    # @todo Should alpha, lambda, and total_effects be logged?
    n_after = (lambdas * density) / (1 + alphas * density + total_effects)

    # Store the updated population densities for this time point in a new data frame.
    # Everything should be in the same order so can index into current density file.
    next_density = pd.DataFrame(columns=density_columns)
    for i in range(len(n_after)):
        if n_after[i] < NO_GROWTH or str(n_after[i]) == 'inf' or str(n_after[i]) == 'nan':
            n_after[i] = 0.0
        next_density = next_density.append({'MODEL_ID': current_density['MODEL_ID'][i], 'DENSITY': n_after[i]},
                                           ignore_index=True)
    next_density.to_csv(density_file_name, index=False)

    return next_density


def get_exchange_fluxes(single_file_names, current_diet, single_rate_file_name, time_point_id,
                        n_processes=None, solver=None):
    """ Optimize the single species models to get exchange reaction fluxes.

    Need more explanation here on why this is done.

    Parameters
    ----------
    single_file_names : list of str
        Paths to single species model files of community members
    current_diet : dict
        Dictionary with exchange reaction ID as key and bound as value
    single_rate_file_name : str
        Path to file for storing growth rates of single species models
    time_point_id : str
        ID of current time point
    n_processes : int, optional
        Number of processes in job pool or None to run without multiprocessing
    solver : str, optional
        Name of solver to use for optimization

    Returns
    -------
    dict
        Dictionary keyed by model ID of fluxes for exchange reactions
    """

    logger.info('[%s] Optimizing %d single species models for exchange fluxes',
                time_point_id, len(single_file_names))

    # Optimize all of the single species models on the current diet conditions.
    if n_processes is None:
        detail_list = [optimize_single_model(file_name, current_diet, solver=solver)
                       for file_name in single_file_names]
    else:
        pool = Pool(n_processes)
        result_list = [pool.apply_async(optimize_single_model, (file_name, current_diet, 'e', solver))
                       for file_name in single_file_names]
        detail_list = [result.get() for result in result_list]
        pool.close()

    # Get the fluxes for the metabolites that are consumed and produced by each species.
    single_rate = pd.DataFrame(columns=single_rate_columns)
    exchange_fluxes = dict()
    for details in detail_list:
        if details['objective_value'] < NO_GROWTH:
            logger.warn('[%s] Single model %s did not grow on current diet conditions: %s %f',
                        time_point_id, details['model_id'], details['status'], details['objective_value'])
        exchange_fluxes[details['model_id']] = details['exchange_fluxes']
        rate = pd.Series([details['model_id'], details['status'], details['objective_value']],
                         index=single_rate_columns)
        single_rate = single_rate.append(rate, ignore_index=True)
    single_rate.to_csv(single_rate_file_name, index=False)
    return exchange_fluxes


def create_next_diet(current_diet, exchange_fluxes, density, next_diet_file_name, time_step, time_point_id):
    """ Create diet conditions for next time point from current diet conditions.
    
    Parameters
    ----------
    current_diet : dict
        Dictionary with exchange reaction ID as key and bound as value
    exchange_fluxes : dict
        Dictionary keyed by model ID of fluxes for exchange reactions
    density : pandas.DataFrame
        Population densities for species in community
    next_diet_file_name: str
        Path to file for storing diet conditions of next time point
    time_step : float
        Need a description here
    time_point_id : str
        ID of current time point

    Returns
    -------
    dict
        Dictionary with exchange reaction ID as key and bound as value
    """

    logger.info('[%s] Calculating diet conditions for next time point', time_point_id)

    # Calculate consumption and production of every metabolite in the diet at
    # this time point, adjusted by the population density and time step.
    new_fluxes = defaultdict(float)
    for species_id in exchange_fluxes:
        row = density.loc[density['MODEL_ID'] == species_id]
        for met_id in exchange_fluxes[species_id]:
            value = exchange_fluxes[species_id][met_id] * row.iloc[0]['DENSITY'] * time_step
            new_fluxes[met_id] += value

    # Update the current diet conditions to create the diet conditions for the
    # next time point.
    next_diet = current_diet
    for met_id in new_fluxes:
        next_diet[met_id] += new_fluxes[met_id]
        if next_diet[met_id] < 0.:
            next_diet[met_id] = 0.
    json.dump(next_diet, open(next_diet_file_name, 'w'), indent=4)

    return next_diet
