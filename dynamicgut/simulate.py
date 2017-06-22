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
pair_rate_columns = ['A_ID', 'B_ID', 'A_TOGETHER', 'A_ALONE', 'A_CHANGE', 'B_TOGETHER', 'B_ALONE', 'B_CHANGE']

# Column names for single species growth rate data frame.
single_rate_columns = ['ID', 'STATUS', 'GROWTH_RATE']

# Column names for population density data frame.
density_columns = ['ID', 'DENSITY']

# Minimum objective value to show growth.
NO_GROWTH = 1e-13

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
                   density_file_name, data_folder, n_processes=None, verbose=True):
    """ Run a simulation over a time interval.

    Do we need a time_step parameter?
    
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
    n_processes: int, optional
        Number of processes in job pool
    verbose : bool
        Store intermediate data generated at each time point
    """

    # Do the single models need to produce growth on given diet?

    # @todo Does time step need to a parameter?
    time_step = 0.5

    # Get the initial population density values.
    # @todo Somewhere need to confirm that IDs in density match IDs in single models
    density = pd.read_csv(density_file_name, dtype={'ID': str, 'DENSITY': float})
    if not set(density_columns).issubset(density.columns):
        raise ValueError('Required columns {0} not available in population density file {1}'
                         .format(density_columns, density_file_name))

    # Get the initial diet conditions.
    initial_diet = json.load(open(diet_file_name))

    # Confirm that the initial diet includes every exchange reaction from the
    # single species models.
    model_exchanges = get_exchange_reaction_ids(single_file_names)
    initial_exchanges = set(initial_diet.keys())
    if not initial_exchanges.issuperset(model_exchanges):
        raise ValueError('Diet file {0} does not contain all exchange reactions from single species models'
                         .format(diet_file_name))
    # @todo Warning if initial exchanges has more reactions than models?

    # Create a job pool for running optimizations.
    if n_processes is None:
        n_processes = min(cpu_count(), 4)
    pool = Pool(n_processes)

    # Confirm input single species models produce growth using initial diet conditions.
    result_list = [pool.apply_async(optimize_single_model, (file_name, initial_diet))
                   for file_name in single_file_names]
    for result in result_list:
        details = result.get()
        if details['objective_value'] < NO_GROWTH:
            raise ValueError('Model {0} does not grow using initial diet'.format(details['model_id']))
    current_diet = initial_diet

    # Run the simulation over the specified time interval.
    for time_point in time_interval:
        # Start this time point.
        time_point_id = '{0:04d}'.format(time_point + 1)
        LOGGER.info('[%s] STARTED TIME POINT', time_point_id)
        time_point_folder = join(data_folder, 'timepoint-'+time_point_id)
        if not exists(time_point_folder):
            makedirs(time_point_folder)
        pair_rate_file_name = join(time_point_folder, 'pair-rates-{0}.csv'.format(time_point_id))
        gr_matrix_filename = join(time_point_folder, 'rates_matrix-{0:04d}.csv'.format(time_point))
        single_rate_file_name = join(time_point_folder, 'single-rates-{0}.csv').format(time_point_id)
        next_diet_file_name = join(time_point_folder, 'diet-{0}.json'.format(time_point_id))

        # Calculate the growth rates for each two species model under the current diet conditions.
        growth_rates = calculate_growth_rates(pair_file_names, current_diet, pool, pair_rate_file_name, time_point_id)

        # Accumulate single species growth rates into a dict (confirm get the same answer every time.

        # effects = growth_rates.apply(get_effects, axis=1)
        # if data_folder:
        #     effects.to_csv(join(time_point_folder, 'effects-{0}.csv'.format(time_point_id)))

        # Create the growth rate matrix.
        gr_matrix = create_gr_matrix(growth_rates)
        gr_matrix.to_csv(gr_matrix_filename)

        # Create the effects matrix.
        effects_matrix = create_effects_matrix(growth_rates)
        effects_matrix_filename = join(time_point_folder, 'effects_matrix-{0}.csv'.format(time_point_id))
        effects_matrix.to_csv(effects_matrix_filename)

        # Run Leslie-Gower algorithm to calculate new population densities.
        LOGGER.info('[%s] Calculating population densities ...', time_point_id)
        density = leslie_gower(gr_matrix_filename, effects_matrix_filename, density)
        if data_folder:
            density.to_csv(join(time_point_folder, 'density-{0}.csv'.format(time_point_id)))

        # Get the exchange reaction fluxes from optimizing single species models.
        exchange_fluxes = get_exchange_fluxes(single_file_names, current_diet, pool,
                                              single_rate_file_name, time_point_id)

        # Create diet conditions for next time point.
        current_diet = create_next_diet(current_diet, exchange_fluxes, density,
                                        next_diet_file_name, time_step, time_point_id)

    pool.close()

    return


def optimize_pair_model(model_file_name, medium):
    """ Optimize a two species community model.

    This function is used as a target function in a multiprocessing pool.

    Note that the model is read from a file each time so there is no need to revert
    the model after the optimization.

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
    if t_solution.status == 'optimal' and a_solution.status == 'optimal' and b_solution.status == 'optimal':
        a_alone = a_solution.fluxes[a_objective]
        if a_alone == 0.0:
            a_alone = float(1e-25)
        a_together = t_solution.fluxes[a_objective]
        a_percent_change = a_together / a_alone

        b_alone = b_solution.fluxes[b_objective]
        if b_alone == 0.0:
            b_alone = float(1e-25)
        b_together = t_solution.fluxes[b_objective]
        b_percent_change = b_together / b_alone
        details = pd.Series([a_id, b_id, a_together, a_alone, a_percent_change,
                             b_together, b_alone, b_percent_change],
                            index=pair_rate_columns)
    else:
        details = pd.Series([a_id, b_id, 0., 0., 0., 0., 0., 0.], index=pair_rate_columns)
        if t_solution.status == 'optimal':
            details.set_value('A_TOGETHER', t_solution.fluxes[a_objective])
            details.set_value('B_TOGETHER', t_solution.fluxes[b_objective])
        if a_solution.status == 'optimal':
            details.set_value('A_ALONE', a_solution.fluxes[a_objective])
        if b_solution.status == 'optimal':
            details.set_value('B_ALONE', b_solution.fluxes[b_objective])
        # @todo Need to rework this to get more data. what if 2 of 3 solutions are optimal?

    return details


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

    # for index, row in growth_rates.iterrows():
    return pair_rate


def create_gr_matrix(growth_rates):
    # Matrix has growth rate of species alone on diagonal, species in the presence of other species
    # in the rest of the cells. Column 0 is an organism ID, row 0 is an organism ID

    # But really just need the growth rate of the organism on the diagonal.
    genomeAdata = []
    genomeBdata = []
    grAmixdata = []
    grBmixdata = []
    grAsolodata = []
    grBsolodata = []
    for index, row in growth_rates.iterrows():
        genomeAdata.append(row['A_ID'])  # A_ID
        genomeBdata.append(row['B_ID'])  # B_ID
        grAmixdata.append(row['A_TOGETHER'])  # A_TOGETHER
        grBmixdata.append(row['B_TOGETHER'])  # B_TOGETHER
        grAsolodata.append(row['A_ALONE'])  # A_ALONE
        grBsolodata.append(row['B_ALONE'])  # B_ALONE
    raw_data_df1 = {'Growth_of': genomeAdata, 'In_presence_of': genomeBdata, 'gr': grAmixdata}
    raw_data_df2 = {'Growth_of': genomeBdata, 'In_presence_of': genomeAdata, 'gr': grBmixdata}
    raw_data_df3 = {'Growth_of': genomeAdata, 'In_presence_of': genomeAdata, 'gr': grAsolodata}
    raw_data_df4 = {'Growth_of': genomeBdata, 'In_presence_of': genomeBdata, 'gr': grBsolodata}
    df1 = pd.DataFrame(raw_data_df1, columns=['Growth_of', 'In_presence_of', 'gr'])
    df2 = pd.DataFrame(raw_data_df2, columns=['Growth_of', 'In_presence_of', 'gr'])
    df3 = pd.DataFrame(raw_data_df3, columns=['Growth_of', 'In_presence_of', 'gr'])
    df4 = pd.DataFrame(raw_data_df4, columns=['Growth_of', 'In_presence_of', 'gr'])
    new_df = df1
    new_df = new_df.append(df2, ignore_index=True)
    new_df = new_df.append(df3, ignore_index=True)
    new_df = new_df.append(df4, ignore_index=True)
    new_df = new_df.pivot_table(index = 'Growth_of', columns = 'In_presence_of', values = 'gr')
    return new_df


def create_effects_matrix(effects):
    genomeAdata = []
    genomeBdata = []
    percentChangeA = []
    percentChangeB = []
    for index, row in effects.iterrows():
        genomeAdata.append(row['A_ID'])
        genomeBdata.append(row['B_ID'])
        percentChangeA.append(row['A_CHANGE'])
        percentChangeB.append(row['B_CHANGE'])
    raw_data_df1 = {'Percent_change_in_growth_of': genomeAdata, 'Because_of': genomeBdata, 'change': percentChangeA}
    raw_data_df2 = {'Percent_change_in_growth_of': genomeBdata, 'Because_of': genomeAdata, 'change': percentChangeB}
    df1 = pd.DataFrame(raw_data_df1, columns=['Percent_change_in_growth_of', 'Because_of', 'change'])
    df2 = pd.DataFrame(raw_data_df2, columns=['Percent_change_in_growth_of', 'Because_of', 'change'])
    new_df = df1
    new_df = new_df.append(df2, ignore_index=True)
    new_df = new_df.pivot_table(index='Percent_change_in_growth_of', columns='Because_of', values='change')
    new_df = new_df.replace(np.nan,'1', regex = True)
    return new_df


def leslie_gower(gr_matrix_filename, effects_matrix_filename, density, k=1, time_step=0.5):
    """ Run Leslie-Gower algorithm to update population density of each organism.
    
    Parameters
    ----------
    gr_matrix_filename : str
        Path to 
    effects_matrix_filename : str
        Path to 
    density : pandas.DataFrame
        Current population density for organisms in community
    k : int
        Something important (upper asymptote in numbers) K = r/a where r is difference between birth-rate 
        and death-rate and a is a positive constant paper tried K of 100 and 2000
    time_step : float
        Size of time interval where 1 means one hour
    """

    # Need to understand the input data. Can it be in a better format?

    # And who is Gower?
    
    # Need consistent order to arrays.
    density_data = []
    for index, row in density.iterrows():
        density_data.append(row['DENSITY'])
    initial_density = np.array(density_data, dtype=float)

    # Gotta be careful that organism IDs all match up.

    # Gotta be a better way to do this but column names are different based on organisms in simulation.
    effects_data = []
    with open(effects_matrix_filename, 'r') as handle:
        handle.readline()
        for line in handle:
            fields = line.strip().split(',')
            effects_data.append(fields[1:])
    effects = np.array(effects_data, dtype=float)
    # So "effects" is a matrix with row 0 and column 0 removed from the input file (which started as a DataFrame
    # It is really important that everything is in the right order

    # Actual Effects (21 January 2017). This is the decrease in growth of the focal species due to the others.
    # Previously, we had the total growth of the focal species in the presence of the other, which i don't
    # think made much sense for being in the numerators.
    # Actually, I'm still not completely sure about this.
    effects = 1 - effects

    # Calculate the vector of the total effects of other species on each of our focal species
    sum_effects = np.dot(effects, initial_density) #Ok, this is the right way.

    # Get the information relative to how much biomass is created per species under the present conditions (Bt)
    species_biomasses = extract_biomass(gr_matrix_filename) #remember that column 1 is the speciesIDs and column 2 is biomasses
    # This is a complicated way to get the growth rate of something?

    #get just the biomasses in a vector
    species_ids = [] # This is used to build the output data frame
    Bt = []  # Birth rate at time t? BetaT What about the death rate DeltaT

    # lambdaT = 1 + BetaT - DeltaT equation 2.3

    # why a population density instead of population size?

    for line in species_biomasses:
        species_ids.append(line[0])
        Bt.append(line[1])

    species_ids = species_ids[1:]
    Bt = Bt[1:]
    Bt = np.array(Bt, dtype=float) # This is a numpy array of single species growth rate on current diet

    #reduce the size of the time step
    # what about when time step is 0?
    Bt = Bt * time_step # birth rate

    # 2 species vs predator-prey

    #create a vector of lambdas
    #FIXME Attention: lbds should be equal to 1 + (Bt over the initial population size used for the
    # calculation of Bt, so 1). This way, Bt will also be independent of populations size
    #The initial calculations were using the below and that may be why everything always seemed so strange
    #lbds = 1 + Bt/init

    # I think it is supposed to be
    lbds = 1 + Bt

    #create a vector of alphas
    alphas = (lbds-1) / k

    #create a vector with the values of the population sizes at Nt+1
    Nafter = (lbds * initial_density) / (1 + alphas * initial_density + sum_effects)

    # Create a new densities file.
    density = pd.DataFrame(columns=['ORGANISM', 'DENSITY'])
    for i in range(len(Nafter)):
        if Nafter[i] < 1e-12 or str(Nafter[i]) == 'inf' or str(Nafter[i]) == 'nan':
            Nafter[i] = 0.0
        density = density.append({'ORGANISM': species_ids[i], 'DENSITY': Nafter[i]}, ignore_index=True)

    return density


def extract_biomass(gr_matrix_filename):

    # Gotta be a better way to do this ...
    matrix = []
    with open(gr_matrix_filename, 'r') as handle:
        for line in handle:
            matrix.append(line.strip().split(','))

    #create an array jsut with the speciesIDs
    speciesIDs = []
    for i in range(len(matrix)):
        speciesIDs.append(matrix[i][0])
    speciesIDs[0] = 'SpeciesIDs'

    #create an array with the biomass created for each species in isolation
    col2 = []
    for i in range(len(matrix)):
        for j in range(len(matrix)):
            if i == j:
                col2.append(matrix[i][j])
    col2[0] = 'Biomass'

    #merge the information matching the biomass values to the speciesID
    merged = []
    for i in range(len(speciesIDs)):
        new_line = [speciesIDs[i],col2[i]]
        merged.append(new_line)

    # So merged is a list of lists, first row is a header with 'SpeciesIDs" and 'Biomasses'
    # all remaining rows are the model ID and growth rate of species by itself in current diet
    return merged


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

    # Get the fluxes for the metabolites that are consumed and produced by each organism.
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
        Population densities for organisms in community
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
    for organism_id in exchange_fluxes:
        row = density.loc[density['ORGANISM']==organism_id]
        for rxn_id in exchange_fluxes[organism_id]:
            value = exchange_fluxes[organism_id][rxn_id] * row.iloc[0]['DENSITY'] * time_step
            new_fluxes[rxn_id] += value

    # Update the current diet conditions to create the diet conditions for the
    # next time point.
    next_diet = current_diet
    for rxn_id in new_fluxes:
        next_diet[rxn_id] += new_fluxes[rxn_id]
    json.dump(next_diet, open(next_diet_file_name, 'w'), indent=4)

    return next_diet
