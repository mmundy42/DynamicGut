from os.path import join
from shutil import rmtree

import dynamicgut

class TestSimulation:

    def test_simulation(self, data_folder, model_files, test_folder):
        single_folder = join(test_folder, 'single')
        pair_folder = join(test_folder, 'pair')
        single_models, pair_models = dynamicgut.prepare_simulation([join(data_folder, m) for m in model_files],
                                                                   single_folder, pair_folder)
        assert len(single_models) == 3
        for filename in single_models:
            assert filename.startswith(single_folder)
        assert len(pair_models) == 3
        for filename in pair_models:
            assert filename.startswith(pair_folder)

        diet_filename = join(data_folder, 'initial_unhealthy_diet.json')
        density_filename = join(data_folder, 'initial_density.csv')
        dynamicgut.run_simulation(range(10), single_models, pair_models, diet_filename, density_filename, test_folder)

        # Check the values in timepoint-0009/density-0009.csv
        # What about values in timepoint-0009/effects-0009.csv and rates-0009.csv
        rmtree(test_folder)
