from os import makedirs
from os.path import join, exists
from warnings import warn
from multiprocessing import Pool, cpu_count
import pandas as pd
import numpy as np
import json
import logging

from mminte import calculate_growth_rates, load_model_from_file, single_species_knockout
from mminte import create_interaction_models, get_all_pairs
from mminte.interaction_worker import compute_growth_rates, growth_rate_columns, apply_medium

from .util import check_for_growth, get_exchange_reaction_ids

# Logger for this module
LOGGER = logging.getLogger(__name__)

interaction_columns = ['A_ID', 'B_ID', 'A_TOGETHER', 'A_ALONE', 'A_CHANGE', 'B_TOGETHER', 'B_ALONE', 'B_CHANGE']

# Column names for population density data frame.
density_columns = ['ID', 'DENSITY']

# Column names for single species growth rate data frame.
single_rate_columns = ['ID', 'RATE']

# Minimum objective value to show growth.
NO_GROWTH = 1e-6

# check on objectives in Agora models and which one got picked when building pair model

# For debugging pandas
# set_option('display.width', 1000)


def prepare(single_file_names, pair_model_folder, optimize=False, n_processes=None):
    """ Prepare for simulation by creating two species community models.

    Parameters
    ----------
    single_file_names : list of str
        Paths to input single species model files
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
    LOGGER.info('Started creating pair models')
    pair_file_names = create_interaction_models(get_all_pairs(single_file_names),
                                                pair_model_folder,
                                                n_processes=n_processes)
    LOGGER.info('Finished creating pair models')
    return pair_file_names


def run_simulation(time_interval, single_file_names, pair_file_names, diet_file_name,
                   density_file_name, data_folder, n_processes=None, verbose=True):
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
    n_processes: int, optional
        Number of processes in job pool
    verbose : bool
        Store intermediate data generated at each time point
    """

    # Do the single models need to produce growth on given diet?

    # Get the initial population density values.
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
        time_point_id = '{0:04d}'.format(time_point)
        LOGGER.info('TIME POINT %s', time_point_id)
        time_point_folder = join(data_folder, 'timepoint-'+time_point_id)
        if not exists(time_point_folder):
            makedirs(time_point_folder)

        # load the diet file once and pass around the dict

        # Optimize the single species models to get exchange reaction fluxes. This is used later
        # to adjust the diet conditions.
        LOGGER.info('Optimizing single species models ...')
        single_rate = pd.DataFrame(columns=single_rate_columns)
        exchange_fluxes = dict()
        result_list = [pool.apply_async(optimize_single_model, (file_name, initial_diet))
                       for file_name in single_file_names]
        for result in result_list:
            details = result.get()
            # @todo Check for no growth here?
            exchange_fluxes[details['model_id']] = details['exchange_fluxes']
            rate = pd.Series([details['model_id'], details['objective_value']], index=single_rate_columns)
            single_rate = single_rate.append(rate, ignore_index=True)
        single_rate.to_csv(join(time_point_folder, 'single-rates-{0}.csv').format(time_point_id), index=False)

        # Calculate the growth rates for each two species model under the current diet conditions.
        LOGGER.info('Optimizing pair models ...')
        medium = json.load(open(diet_file_name))
        # Think need to use optimize_pair_model here, and just return the pd.Series that is need for LG calculation
        growth_rates = calculate_growth_rates(pair_file_names, medium)
        if verbose:
            growth_rates.to_csv(join(time_point_folder, 'rates-{0:04d}.csv'.format(time_point)), index=False)

        effects = growth_rates.apply(get_effects, axis=1)
        if data_folder:
            effects.to_csv(join(time_point_folder, 'effects-{0:04d}.csv'.format(time_point)))

        # Create the growth rate matrix.
        gr_matrix = create_gr_matrix(growth_rates)
        gr_matrix_filename = join(time_point_folder, 'rates_matrix-{0:04d}.csv'.format(time_point))
        gr_matrix.to_csv(gr_matrix_filename)

        # Create the effects matrix.
        effects_matrix = create_effects_matrix(effects)
        effects_matrix_filename = join(time_point_folder, 'effects_matrix-{0:04d}.csv'.format(time_point))
        effects_matrix.to_csv(effects_matrix_filename)

        # Run Leslie-Gower algorithm to calculate new population densities.
        LOGGER.info('Calculating population densities ...')
        density = leslie_gower(gr_matrix_filename, effects_matrix_filename, density)
        if data_folder:
            density.to_csv(join(time_point_folder, 'density-{0:04d}.csv'.format(time_point)))

        LOGGER.info('Calculating diet conditions for next time point ...')
        next_diet_filename = join(time_point_folder, 'diet-{0}.json'.format(time_point_id))
        create_next_diet(diet_file_name, next_diet_filename, exchange_fluxes, density)
        diet_file_name = next_diet_filename

    pool.close()

    return


def optimize_single_model(model_file_name, medium):
    """ Optimize a single species model on a given medium.
    
    Note that we chose to read the model from a file each time instead of loading
    the model into memory once at the beginning of the simulation. This lowers
    the memory requirements of the simulation but is slower because of repeated 
    file system accesses.
    
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

    model = load_model_from_file(model_file_name)
    details = {'model_id': model.id}
    apply_medium(model, medium)
    solution = model.optimize()
    details['objective_value'] = solution.objective_value
    exchange_reactions = model.reactions.query(lambda x: x.startswith('EX_'), 'id')
    details['exchange_fluxes'] = dict()
    if solution.status == 'optimal':
        for rxn in exchange_reactions:
            details['exchange_fluxes'][rxn.id] = solution.fluxes[rxn.id]
    else:
        for rxn in exchange_reactions:
            details['exchange_fluxes'][rxn.id] = 0.0
    return details


def optimize_pair_model(pair_filename, media_filename):
    """ Optimize a two species community model.

    Parameters
    ----------
    pair_filename : str
        Path to two species community model file
    media_filename : str
        Path to file with exchange reaction bounds for media

    Returns
    -------
    pandas.Series
        Growth rate details for interaction between two species in pair
    """

    # Load the model and apply the media to it.
    pair_model = load_model_from_file(pair_filename)
    apply_medium(pair_model, media_filename)

    # Optimize the model with two species together, one species knocked out, and
    # other species knocked out.
    result = dict()
    result['a_id'] = pair_model.notes['species'][0]['id']
    result['a_objective'] = pair_model.notes['species'][0]['objective']
    result['b_id'] = pair_model.notes['species'][1]['id']
    result['b_objective'] = pair_model.notes['species'][1]['objective']
    result['t_solution'] = pair_model.optimize()
    result['a_solution'] = single_species_knockout(pair_model, result['b_id'])
    result['b_solution'] = single_species_knockout(pair_model, result['a_id'])

    return result


def get_effects(row):
    """
    """

    a_alone = row['A_ALONE']
    if a_alone == 0.0:
        a_alone = float(1e-25)
    a_percent_change = row['A_TOGETHER'] / a_alone
    b_alone = row['B_ALONE']
    if b_alone == 0.0:
        b_alone = float(1e-25)
    b_percent_change = row['B_TOGETHER'] / b_alone
    details = pd.Series([row['A_ID'], row['B_ID'], row['A_TOGETHER'], row['A_ALONE'], a_percent_change,
                         row['B_TOGETHER'], row['B_ALONE'], b_percent_change], index=interaction_columns)
    return details

    # if results['t_solution'].status == 'optimal' and \
    #     results['a_solution'].status == 'optimal' and \
    #     results['b_solution'].status == 'optimal':
    #     a_alone = results['a_solution'].x_dict[results['a_objective']]
    #     if a_alone == 0.0:
    #         a_alone = float(1e-25)
    #     a_together = results['t_solution'].x_dict[results['a_objective']]
    #     a_percent_change = a_together / a_alone
    #
    #     b_alone = results['b_solution'].x_dict[results['b_objective']]
    #     if b_alone == 0.0:
    #         b_alone = float(1e-25)
    #     b_together = results['t_solution'].x_dict[results['b_objective']]
    #     b_percent_change = b_together / b_alone
    #     details = Series([results['a_id'], results['b_id'], a_together, a_alone, a_percent_change,
    #                       b_together, b_alone, b_percent_change], index=interaction_columns)
    # else:
    #     details = Series([results['a_id'], results['b_id'], 'None', 0., 0., 0., 0., 0., 0., 0.],
    #                      index=growth_rate_columns)
    #     if results['t_solution'].status == 'optimal':
    #         details.set_value('TOGETHER', results['t_solution'].f)
    #         details.set_value('A_TOGETHER', results['t_solution'].x_dict[results['a_objective']])
    #         details.set_value('B_TOGETHER', results['t_solution'].x_dict[results['b_objective']])
    #     if results['a_solution'].status == 'optimal':
    #         details.set_value('A_ALONE', results['a_solution'].x_dict[results['a_objective']])
    #     if results['b_solution'].status == 'optimal':
    #         details.set_value('B_ALONE', results['b_solution'].x_dict[results['b_objective']])

    return details


def create_gr_matrix(growth_rates):
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
        Something important
    time_step : float
        Size of time interval where 1 means one hour
    """

    # Need to understand the input data. Can it be in a better format?
    
    # Need consistent order to arrays.
    density_data = []
    for index, row in density.iterrows():
        density_data.append(row['DENSITY'])
    initial_density = np.array(density_data, dtype=float)

    # Gotta be a better way to do this but column names are different based on organisms in simulation.
    effects_data = []
    with open(effects_matrix_filename, 'r') as handle:
        handle.readline()
        for line in handle:
            fields = line.strip().split(',')
            effects_data.append(fields[1:])
    effects = np.array(effects_data, dtype=float)

    # Actual Effects (21 January 2017). This is the decrease in growth of the focal species due to the others.
    # Previously, we had the total growth of the focal species in the presence of the other, which i don't
    # think made much sense for being in the numerators.
    # Actually, I'm still not completely sure about this.
    effects = 1 - effects

    # Calculate the vector of the total effects of other species on each of our focal species
    sum_effects = np.dot(effects, initial_density) #Ok, this is the right way.

    # Get the information relative to how much biomass is created per species under the present conditions (Bt)
    species_biomasses = extract_biomass(gr_matrix_filename) #remember that column 1 is the speciesIDs and column 2 is biomasses

    #get just the biomasses in a vector
    species_ids = []
    Bt = []

    for line in species_biomasses:
        species_ids.append(line[0])
        Bt.append(line[1])

    species_ids = species_ids[1:]
    Bt = Bt[1:]
    Bt = np.array(Bt, dtype=float)

    #reduce the size of the time step
    # what about when time step is 0?
    Bt = Bt * time_step # birth rate

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

    return merged


def create_next_diet(current_diet_filename, next_diet_filename, exchange_fluxes, density, time_step=0.5):
    """ Create a diet file by updating the current diet.
    
    Parameters
    ----------
    current_diet_filename : str
        X
    next_diet_filename : str
        Path to file
    exchange_fluxes : dict
        Dictionary keyed by organism ID with exchange fluxes
    density : pandas.DataFrame
        Population densities
    time_step : float
        X
    """

    next_medium = json.load(open(current_diet_filename))

    # Really need to pay attention to signs here

    # new dict with total exchange reaction fluxes
    # go through each organism exchange reaction flux
    # for the flux multiple by organism density and time step
    new_fluxes = dict()
    for organism_id in exchange_fluxes:
        row = density.loc[density['ORGANISM']==organism_id]

        for reaction_id in exchange_fluxes[organism_id]:
            value = exchange_fluxes[organism_id][reaction_id] * row.iloc[0]['DENSITY'] * time_step
            try:
                new_fluxes[reaction_id] += value
            except KeyError:
                new_fluxes[reaction_id] = value

    for reaction_id in new_fluxes:
        value = next_medium[reaction_id] - new_fluxes[reaction_id]
        next_medium[reaction_id] = value

    json.dump(next_medium, open(next_diet_filename, 'w'))

    return


# def saved():
#     for pair_filename in pair_models:
#         result = optimize_pair_model(pair_filename, diet_file)
#
#         # Go through all of the exchange fluxes in the two single species solutions.
#         if result['a_id'] not in exchange_fluxes:
#             exchange_fluxes[result['a_id']] = dict()
#         # Will need to confirm every single species knockout has same solver status.
#         if result['a_solution'].status == 'optimal':
#             for reaction_id in result['a_solution'].x_dict:
#                 if reaction_id.startswith('EX_'):
#                     if reaction_id in exchange_fluxes[result['a_id']]:
#                         if result['a_solution'].x_dict[reaction_id] != exchange_fluxes[result['a_id']][reaction_id]:
#                             stop = 1
#                             # warn('{0} {1} != {2}'.format(result['a_id'], result['a_solution'].x_dict[reaction_id],
#                             #      exchange_fluxes[result['a_id']][reaction_id]))
#                     else:
#                         exchange_fluxes[result['a_id']][reaction_id] = result['a_solution'].x_dict[reaction_id]
#         if result['b_id'] not in exchange_fluxes:
#             exchange_fluxes[result['b_id']] = dict()
#         if result['b_solution'].status == 'optimal':
#             for reaction_id in result['b_solution'].x_dict:
#                 if reaction_id.startswith('EX_'):
#                     if reaction_id in exchange_fluxes[result['b_id']]:
#                         if result['b_solution'].x_dict[reaction_id] != exchange_fluxes[result['b_id']][reaction_id]:
#                             stop = 1
#                             # warn('{0} {1} != {2}'.format(result['b_id'],
#                             #                              result['b_solution'].x_dict[reaction_id],
#                             #      exchange_fluxes[result['b_id']][reaction_id]))
#                     else:
#                         exchange_fluxes[result['b_id']][reaction_id] = result['b_solution'].x_dict[reaction_id]
#
#         effects = effects.append(get_effects(result), ignore_index=True)
#
#         # This can be simplified but no time now
#         if result['t_solution'].status == 'optimal' and \
#                         result['a_solution'].status == 'optimal' and \
#                         result['b_solution'].status == 'optimal':
#
#             rates = Series([result['a_id'], result['b_id'], 'unknown', result['t_solution'].f,
#                             result['t_solution'].x_dict[result['a_objective']],
#                             result['t_solution'].x_dict[result['b_objective']],
#                             result['a_solution'].x_dict[result['a_objective']],
#                             result['b_solution'].x_dict[result['b_objective']],
#                             0., 0.], index=growth_rate_columns)
#         else:
#             rates = Series([result['a_id'], result['b_id'], 'unknown', 0., 0., 0., 0., 0., 0., 0.],
#                            index=growth_rate_columns)
#             if result['t_solution'].status == 'optimal':
#                 rates.set_value('TOGETHER', result['t_solution'].f)
#                 rates.set_value('A_TOGETHER', result['t_solution'].x_dict[result['a_objective']])
#                 rates.set_value('B_TOGETHER', result['t_solution'].x_dict[result['b_objective']])
#             if result['a_solution'].status == 'optimal':
#                 rates.set_value('A_ALONE', result['a_solution'].x_dict[result['a_objective']])
#             if result['b_solution'].status == 'optimal':
#                 rates.set_value('B_ALONE', result['b_solution'].x_dict[result['b_objective']])
#         growth_rates = growth_rates.append(rates, ignore_index=True)
#
