from os import makedirs
from os.path import join, exists
from warnings import warn
from multiprocessing import Pool, cpu_count
from pandas import DataFrame, Series, set_option, read_csv
import numpy as np
import json

from mminte import calculate_growth_rates, load_model_from_file, single_species_knockout
from mminte.interaction_worker import compute_growth_rates, growth_rate_columns, apply_medium

interaction_columns = ['A_ID', 'B_ID', 'A_TOGETHER', 'A_ALONE', 'A_CHANGE', 'B_TOGETHER', 'B_ALONE', 'B_CHANGE']

# Column names for population density data frame.
density_columns = ['ID', 'DENSITY']

def x_test(model, medium_file):
    apply_media(model, medium_file)

# check on objectives in Agora models and which one got picked when building pair model

# Density file format:
# CSV with two columns: ID, DENSITY

# Diet file format:
# JSON with exchange reaction ID as key, initial concentration as value

# For debugging pandas
# set_option('display.width', 1000)


def run_simulation(time_interval, single_file_names, pair_file_names, diet_file_name,
                   density_file_name, data_folder, verbose=True):
    """ Run a simulation over a time interval.

    Parameters
    ----------
    time_interval : range
        Range of time points for running the simulation
    single_file_names : list of str
        List of path names to single species models with all exchange reactions
    pair_file_names : list of str
        List of path names to two species community model files
    diet_file_name : str
        Path to file with initial diet conditions in JSON format
    density_file_name : str
        Path to file with initial population densities in CSV format
    data_folder : str
        Path to folder for storing data generated at each time point
    verbose : bool
        Store intermediate data generated at each time point
    """

    # Get the initial population density values.
    density = read_csv(density_file_name)
    if not set(density_columns).issubset(density.columns):
        raise ValueError('Required columns {0} not available in population density file {1}'
                         .format(density_columns, density_file_name))

    # Load all of the single species models into memory.
    # @todo Will this consume too much memory for a large community?
    single_models = [load_model_from_file(model) for model in single_file_names]

    # @todo Should there be a check to make sure single species models have exact match
    #   of exchange reactions to diet file?

    # Probably should pull apply_medium to here so easier to switch to new cobrapy method.

    # Run the simulation over the specified time interval.
    for time_point in time_interval:
        # Start this time point.
        # print('TIME POINT {0:04d}\n---------------'.format(time_point))
        time_point_folder = join(data_folder, 'timepoint-{0:04d}'.format(time_point))
        if not exists(time_point_folder):
            makedirs(time_point_folder)

        # Optimize the single species models to get exchange reaction fluxes. This is used later
        # to adjust the diet conditions.
        exchange_fluxes = dict()
        # so apply the diet, then optimize, go through all of the exchange fluxes, build dict, and put in exchange_fluxes
        for model in single_models:
            exchange_fluxes[model.id] = optimize_single_model(model, diet_file)

        # Calculate the growth rates for each two species model under the current diet conditions.
        growth_rates = calculate_growth_rates(pair_file_names, diet_file)
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
        density = leslie_gower(gr_matrix_filename, effects_matrix_filename, density)
        if data_folder:
            density.to_csv(join(time_point_folder, 'density-{0:04d}.csv'.format(time_point)))

        next_diet_filename = join(time_point_folder, 'diet-{0:04d}.json'.format(time_point))
        create_next_diet(diet_file, next_diet_filename, exchange_fluxes, density)
        diet_file = next_diet_filename
    return



def optimize_single_model(model, medium_filename):

    # Optimize the single species model on the medium.
    apply_medium(model, medium_filename)
    solution = model.optimize()

    # Get the fluxes for the exchange reactions.
    exchange_fluxes = dict()
    if solution.status == 'optimal':
        for reaction_id in solution.x_dict:
            if reaction_id.startswith('EX_'):
                exchange_fluxes[reaction_id] = solution.x_dict[reaction_id]
    else:
        exchange_reactions = model.reactions.query(lambda x: x.startswith('EX_'))
        for r in exchange_reactions:
            exchange_fluxes[r.id] = 0.0
        # make them zero
    return exchange_fluxes


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
    details = Series([row['A_ID'], row['B_ID'], row['A_TOGETHER'], row['A_ALONE'], a_percent_change,
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
    df1 = DataFrame(raw_data_df1, columns=['Growth_of', 'In_presence_of', 'gr'])
    df2 = DataFrame(raw_data_df2, columns=['Growth_of', 'In_presence_of', 'gr'])
    df3 = DataFrame(raw_data_df3, columns=['Growth_of', 'In_presence_of', 'gr'])
    df4 = DataFrame(raw_data_df4, columns=['Growth_of', 'In_presence_of', 'gr'])
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
    df1 = DataFrame(raw_data_df1, columns=['Percent_change_in_growth_of', 'Because_of', 'change'])
    df2 = DataFrame(raw_data_df2, columns=['Percent_change_in_growth_of', 'Because_of', 'change'])
    new_df = df1
    new_df = new_df.append(df2, ignore_index=True)
    new_df = new_df.pivot_table(index='Percent_change_in_growth_of', columns='Because_of', values='change')
    new_df = new_df.replace(np.nan,'1', regex = True)
    return new_df


def leslie_gower(gr_matrix_filename, effects_matrix_filename, density, k = 1, time_step=0.5):
    # time step allows you to define size of time interval (1 hour, 30 minutes, etc)

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
    density = DataFrame(columns=['ORGANISM', 'DENSITY'])
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

    next_medium = json.load(open(current_diet_filename))

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
