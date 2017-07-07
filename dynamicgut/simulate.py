from os import makedirs
from os.path import join, exists
from warnings import warn
from multiprocessing import Pool, cpu_count
import pandas as pd
import numpy as np
import json
import logging
from collections import defaultdict

from mminte import load_model_from_file, single_species_knockout
from mminte import create_interaction_models, get_all_pairs
from mminte.interaction_worker import compute_growth_rates, growth_rate_columns, apply_medium

from .util import check_for_growth, get_exchange_reaction_ids

# Logger for this module
LOGGER = logging.getLogger(__name__)

# Column names for two species community growth rate data frame.
pair_rate_columns = ['A_ID', 'B_ID', 'A_TOGETHER', 'A_ALONE', 'A_EFFECT', 'B_TOGETHER', 'B_ALONE', 'B_EFFECT']

# Column names for single species growth rate data frame.
single_rate_columns = ['ID', 'STATUS', 'GROWTH_RATE']

# Column names for population density data frame.
density_columns = ['ID', 'DENSITY']

# Minimum objective value to show growth.
NO_GROWTH = 1e-13

# Very small number to prevent division by zero.
ALMOST_ZERO = 1e-25

# check on objectives in Agora models and which one got picked when building pair model

# For debugging pandas
# set_option('display.width', 1000)


def prepare(single_file_names, pair_model_folder, optimize=False, n_processes=None):
    """ Prepare for simulation by creating two species community models.

    Parameters
    ----------
    single_file_names : list of str
        Paths to single species model files of community members
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

    # If requested, confirm all of the single species input models produce growth as provided.
    if optimize:
        LOGGER.info('Started checking %d single models for growth', len(single_file_names))
        if n_processes is None:
            n_processes = min(cpu_count(), 4)
        pool = Pool(n_processes)
        result_list = [pool.apply_async(check_for_growth, (file_name,))
                       for file_name in single_file_names]
        for result in result_list:
            summary = result.get()
            if not summary['grows']:
                warn(summary['message'])
        LOGGER.info('Finished checking single models for growth')

    # Create folder for output pair models if needed.
    if not exists(pair_model_folder):
        makedirs(pair_model_folder)

    # Create all of the pair models and store in specified folder.
    LOGGER.info('Started creating two species community models')
    pair_file_names = create_interaction_models(get_all_pairs(single_file_names),
                                                pair_model_folder,
                                                n_processes=n_processes)
    LOGGER.info('Finished creating two species community models')
    return pair_file_names


def run_simulation(time_interval, single_file_names, pair_file_names, diet_file_name,
                   density_file_name, data_folder, time_step=0.5, k=1, n_processes=None):
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
        Number of processes in job pool
    """

    # Validate time_step parameter.
    if time_step <= 0.0 or time_step > 1.0:
        raise ValueError('time_step parameter must be a value greater than 0 and less than or equal to 1')

    # Get the initial population density values.
    density = pd.read_csv(density_file_name, dtype={'ID': str, 'DENSITY': float})
    if not set(density_columns).issubset(density.columns):
        raise ValueError('Required columns {0} not available in population density file {1}'
                         .format(density_columns, density_file_name))
    # @todo Add more validation of density file, ID is string greater than 0, DENSITY is valid float

    # Get the initial diet conditions.
    diet = json.load(open(diet_file_name))

    # Set diet for first time step by adding exchange reactions from the single
    # species models that are not in the initial diet. This allows metabolites
    # produced by a species to become available in the diet conditions during the
    # simulation.
    model_exchanges, model_ids = get_exchange_reaction_ids(single_file_names)
    initial_exchanges = set(diet.keys())
    if initial_exchanges > model_exchanges:
        warn('Diet file {0} contains more exchange reactions than there are in single species models'
             .format(diet_file_name))
    if model_exchanges.issuperset(initial_exchanges):  # @todo is this necessary?
        for rxn_id in (model_exchanges - initial_exchanges):
            diet[rxn_id] = 0.0
    json.dump(diet, open(join(data_folder, 'initial-diet.json'), 'w'), indent=4)

    # Confirm the model IDs in the initial density file match the model IDs in the
    # list of single species models.
    if density.shape[0] != len(model_ids):
        raise ValueError('Number of species ({0}) in initial density file does not match '
                         'number of single species models ({1})'.format(density.shape[0], len(model_ids)))
    if density.loc[density['ID'].isin(model_ids)].shape[0] != len(model_ids):
        # @todo Figure what does not match and report to caller
        raise ValueError('One or more model IDs in initial density file do not match '
                         'model IDs in single species models')

    # Create a job pool for running optimizations.
    if n_processes is None:
        n_processes = min(cpu_count(), 4)
    pool = Pool(n_processes)

    # Run the simulation over the specified time interval.
    for time_point in time_interval:
        # Start this time point.
        time_point_id = '{0:04d}'.format(time_point + 1)
        LOGGER.info('[%s] STARTED TIME POINT', time_point_id)
        time_point_folder = join(data_folder, 'timepoint-'+time_point_id)
        if not exists(time_point_folder):
            makedirs(time_point_folder)
        pair_rate_file_name = join(time_point_folder, 'pair-rates-{0}.csv'.format(time_point_id))
        effects_matrix_file_name = join(time_point_folder, 'effects-matrix-{0}.csv'.format(time_point_id))
        density_file_name = join(time_point_folder, 'density-{0}.csv'.format(time_point_id))
        single_rate_file_name = join(time_point_folder, 'single-rates-{0}.csv').format(time_point_id)
        next_diet_file_name = join(time_point_folder, 'diet-{0}.json'.format(time_point_id))

        # Calculate the growth rates for each two species model under the current diet conditions.
        growth_rates, alone = calculate_growth_rates(pair_file_names, diet, pool, pair_rate_file_name, time_point_id)

        # Create the effects matrix.
        effects_matrix = create_effects_matrix(growth_rates, effects_matrix_file_name, time_point_id)

        # Run Leslie-Gower algorithm to calculate new population densities.
        density = leslie_gower(effects_matrix, density, density_file_name, time_point_id, alone, k, time_step)

        # Get the exchange reaction fluxes from optimizing single species models.
        exchange_fluxes = get_exchange_fluxes(single_file_names, diet, pool,
                                              single_rate_file_name, time_point_id)

        # Create diet conditions for next time point.
        diet = create_next_diet(diet, exchange_fluxes, density, next_diet_file_name, time_step, time_point_id)

    # Cleanup and store results from last time step in data folder.
    pool.close()
    json.dump(diet, open(join(data_folder, 'final-diet.json'), 'w'), indent=4)
    density.to_csv(join(data_folder, 'final-density.csv'))

    return


def optimize_pair_model(model_file_name, medium):
    """ Optimize a two species community model.

    This function is used as a target function in a multiprocessing pool. Since the
    model is read from a file each time the function runs there is no need to revert
    the model after the optimization.

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
        Dictionary with exchange reaction ID as key and bound as value

    Returns
    -------
    pandas.Series
        Growth rate details for interaction between two species in pair
    """

    # Optimize the model with two species together, one species knocked out, and
    # other species knocked out.
    pair_model = load_model_from_file(model_file_name)
    apply_medium(pair_model, medium)

    a_id = pair_model.notes['species'][0]['id']
    a_objective = pair_model.notes['species'][0]['objective']
    b_id = pair_model.notes['species'][1]['id']
    b_objective = pair_model.notes['species'][1]['objective']

    t_solution = pair_model.optimize()
    a_solution = single_species_knockout(pair_model, b_id)
    b_solution = single_species_knockout(pair_model, a_id)

    # Round very small growth rates to zero.
    if t_solution.fluxes[a_objective] < NO_GROWTH:
        t_solution.fluxes[a_objective] = 0.
    if t_solution.fluxes[b_objective] < NO_GROWTH:
        t_solution.fluxes[b_objective] = 0.
    if a_solution.fluxes[a_objective] < NO_GROWTH:
        a_solution.fluxes[a_objective] = 0.
    if b_solution.fluxes[b_objective] < NO_GROWTH:
        b_solution.fluxes[b_objective] = 0.

    # Evaluate the interaction between the two species.
    if t_solution.status == 'optimal':
        a_together = t_solution.fluxes[a_objective]
        b_together = t_solution.fluxes[b_objective]
    else:
        a_together = 0.0
        b_together = 0.0

    if a_solution.status == 'optimal':
        a_alone = a_solution.fluxes[a_objective]
    else:
        a_alone = 0.0
    if a_alone != 0.0:
        alone = a_alone
    else:
        alone = ALMOST_ZERO
    if a_together != 0.0 or a_alone != 0.0:
        a_effect = a_together / alone  # See note above for description
    else:
        a_effect = 0.0

    if b_solution.status == 'optimal':
        b_alone = b_solution.fluxes[b_objective]
    else:
        b_alone = 0.0
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


def calculate_growth_rates(pair_file_names, current_diet, pool, pair_rate_file_name, time_point_id):
    """ Optimize the two species community models to get growth rates.

    Parameters
    ----------
    pair_file_names : list of str
        Paths to two species community model files
    current_diet : dict
        Dictionary with exchange reaction ID as key and bound as value
    pool : multiprocessing.Pool
        Job pool for running optimizations
    pair_rate_file_name : str
        Path to file for storing growth rates of two species community models
    time_point_id : str
        ID of current time point

    Returns
    -------
    pandas.DataFrame
        Growth rates
    dict
        Single species growth rates
    """

    LOGGER.info('[%s] Optimizing two species community models ...', time_point_id)

    # Optimize all of the two species community models on the current diet conditions.
    result_list = [pool.apply_async(optimize_pair_model, (file_name, current_diet))
                   for file_name in pair_file_names]

    # Build a DataFrame with the pair growth rates.
    pair_rate = pd.DataFrame(columns=pair_rate_columns)
    for result in result_list:
        pair_rate = pair_rate.append(result.get(), ignore_index=True)
    pair_rate.to_csv(pair_rate_file_name, index=False)

    # Build a dictionary with single species growth rates and confirm that the
    # values are consistent.
    alone_rate = dict()
    for index, row in pair_rate.iterrows():
        if row['A_ID'] in alone_rate:
            if not np.isclose(row['A_ALONE'], alone_rate[row['A_ID']]):
                warn('Model {0} has inconsistent growth rates: {1} vs {2}'
                     .format(row['A_ID'], row['A_ALONE'], alone_rate[row['A_ID']]))
        else:
            alone_rate[row['A_ID']] = row['A_ALONE']
        if row['B_ID'] in alone_rate:
            if not np.isclose(row['B_ALONE'], alone_rate[row['B_ID']]):
                warn('Model {0} has inconsistent growth rates: {1} vs {2}'
                     .format(row['B_ID'], row['B_ALONE'], alone_rate[row['B_ID']]))
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

    LOGGER.info('[%s] Creating effects matrix ...', time_point_id)

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

    LOGGER.info('[%s] Calculating population densities ...', time_point_id)

    # Confirm that the order of the species in the effects matrix matches the
    # order in the density data frame.
    if not np.array_equal(effects_matrix.axes[0].values, current_density['ID'].values):
        # @todo Figure out a way to include the mismatched data in the exception or log it
        raise ValueError('Order of species in effects matrix does not match density data frame')

    # The "N1(t)" term in equation above is the population size at the beginning of the
    # time point and comes from the DENSITY column of the input data frame. Create a
    # 1-D array from the data frame.
    density = current_density['DENSITY'].values

    # Convert the dictionary of single species growth rates at this time point to
    # a 1-D array in the same order as the current density array.
    growth = np.array([alone_rate[row['ID']] for index, row in current_density.iterrows()], dtype=float)

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
    effects = 1 - effects

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
        next_density = next_density.append({'ID': current_density['ID'][i], 'DENSITY': n_after[i]}, ignore_index=True)
    next_density.to_csv(density_file_name, index=False)

    return next_density


def optimize_single_model(model_file_name, medium):
    """ Optimize a single species model on a given medium.

    This function is used as a target function in a multiprocessing pool.

    Note that we chose to read the model from a file each time instead of loading
    the model into memory once at the beginning of the simulation. This lowers
    the memory requirements of the simulation and there is no need to revert the
    model after the optimization. But there are more accesses of the file system.

    Parameters
    ----------
    model_file_name : cobra.core.Model
        Single species model to be optimized
    medium : dict
        Dictionary with exchange reaction ID as key and bound as value

    Returns
    -------
    dict
        Dictionary with details on solution
    """

    # Confirmed that growth rates are the same as when run solo in pair model.
    # Are we only doing this to get the exchange reaction fluxes which are
    # unavailable from mminte output?

    model = load_model_from_file(model_file_name)
    details = {'model_id': model.id}
    apply_medium(model, medium)
    solution = model.optimize()
    details['status'] = solution.status
    details['objective_value'] = solution.objective_value
    exchange_reactions = model.reactions.query(lambda x: x.startswith('EX_'), 'id')
    details['exchange_fluxes'] = dict()
    if solution.status == 'optimal':
        for rxn in exchange_reactions:
            if solution.fluxes[rxn.id] != 0.0:
                details['exchange_fluxes'][rxn.id] = solution.fluxes[rxn.id]
    return details


def get_exchange_fluxes(single_file_names, current_diet, pool, single_rate_file_name, time_point_id):
    """ Optimize the single species models to get exchange reaction fluxes.

    Need more explanation here on why this is done.

    Parameters
    ----------
    single_file_names : list of str
        Paths to single species model files of community members
    current_diet : dict
        Dictionary with exchange reaction ID as key and bound as value
    pool : multiprocessing.Pool
        Job pool for running optimizations
    single_rate_file_name : str
        Path to file for storing growth rates of single species models
    time_point_id : str
        ID of current time point

    Returns
    -------
    dict
        Dictionary keyed by model ID of fluxes for exchange reactions
    """

    LOGGER.info('[%s] Optimizing single species models ...', time_point_id)

    # Optimize all of the single species models on the current diet conditions.
    single_rate = pd.DataFrame(columns=single_rate_columns)
    result_list = [pool.apply_async(optimize_single_model, (file_name, current_diet))
                   for file_name in single_file_names]

    # Get the fluxes for the metabolites that are consumed and produced by each species.
    exchange_fluxes = dict()
    for result in result_list:
        details = result.get()
        if details['objective_value'] < NO_GROWTH:
            LOGGER.warn('[%s] Model %s did not grow on current diet conditions', time_point_id, details['model_id'])
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

    LOGGER.info('[%s] Calculating diet conditions for next time point ...', time_point_id)

    # Calculate consumption and production of every metabolite in the diet at
    # this time point, adjusted by the population density and time step.
    new_fluxes = defaultdict(float)
    for species_id in exchange_fluxes:
        row = density.loc[density['ID'] == species_id]
        for rxn_id in exchange_fluxes[species_id]:
            value = exchange_fluxes[species_id][rxn_id] * row.iloc[0]['DENSITY'] * time_step
            new_fluxes[rxn_id] += value

    # Update the current diet conditions to create the diet conditions for the
    # next time point.
    next_diet = current_diet
    for rxn_id in new_fluxes:
        next_diet[rxn_id] += new_fluxes[rxn_id]
    json.dump(next_diet, open(next_diet_file_name, 'w'), indent=4)

    return next_diet
