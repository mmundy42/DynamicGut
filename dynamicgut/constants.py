# Column names for single species growth rate data frame.
single_rate_columns = ['MODEL_ID', 'STATUS', 'GROWTH_RATE']

# Column names for population density data frame.
density_columns = ['MODEL_ID', 'DENSITY']

# Column names for two species community growth rate data frame.
pair_rate_columns = ['A_ID', 'B_ID', 'STATUS', 'A_TOGETHER', 'A_ALONE', 'A_EFFECT', 'B_TOGETHER', 'B_ALONE', 'B_EFFECT']

# Minimum objective value to show growth.
NO_GROWTH = 1e-13

# Very small number to prevent division by zero.
ALMOST_ZERO = 1e-25

# Number of seconds to allow solver to run
SOLVER_TIME_LIMIT = 30
