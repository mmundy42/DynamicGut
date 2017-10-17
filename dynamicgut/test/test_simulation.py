from os.path import join, exists
from shutil import rmtree
import pytest

import dynamicgut


class TestSimulation:

    def test_simulation(self, data_folder, model_files, pair_model_files, test_folder):
        single_models = [join(data_folder, m) for m in model_files]
        pair_models = [join(data_folder, m) for m in pair_model_files]
        diet_file = join(data_folder, 'initial_unhealthy_diet.json')
        density_file = join(data_folder, 'initial_density.csv')
        dynamicgut.run_simulation(range(3), single_models, pair_models, diet_file, density_file, test_folder)
        for index in range(1, 4):
            timepoint = 'timepoint-{0:04d}'.format(index)
            assert exists(join(test_folder, timepoint, 'density-{0:04d}.csv'.format(index)))
            assert exists(join(test_folder, timepoint, 'diet-{0:04d}.json'.format(index)))
            assert exists(join(test_folder, timepoint, 'effects-matrix-{0:04d}.csv'.format(index)))
            assert exists(join(test_folder, timepoint, 'pair-rates-{0:04d}.csv'.format(index)))
            assert exists(join(test_folder, timepoint, 'single-rates-{0:04d}.csv'.format(index)))
        rmtree(test_folder)

    def test_bad_single_path(self, data_folder, model_files, pair_model_files, test_folder):
        single_models = [join(data_folder, m) for m in model_files]
        single_models.append(join(data_folder, 'bad.xml'))
        pair_models = [join(data_folder, m) for m in pair_model_files]
        diet_file = join(data_folder, 'initial_unhealthy_diet.json')
        density_file = join(data_folder, 'initial_density.csv')
        with pytest.raises(IOError):
            dynamicgut.run_simulation(range(3), single_models, pair_models, diet_file, density_file, test_folder)

    def test_bad_pair_path(self, data_folder, model_files, pair_model_files, test_folder):
        single_models = [join(data_folder, m) for m in model_files]
        pair_models = [join(data_folder, m) for m in pair_model_files]
        pair_models.append(join(data_folder, 'bad.pickle'))
        diet_file = join(data_folder, 'initial_unhealthy_diet.json')
        density_file = join(data_folder, 'initial_density.csv')
        with pytest.raises(IOError):
            dynamicgut.run_simulation(range(3), single_models, pair_models, diet_file, density_file, test_folder)

    def test_bad_diet_path(self, data_folder, model_files, pair_model_files, test_folder):
        single_models = [join(data_folder, m) for m in model_files]
        pair_models = [join(data_folder, m) for m in pair_model_files]
        diet_file = join(data_folder, 'bad.json')
        density_file = join(data_folder, 'initial_density.csv')
        with pytest.raises(IOError):
            dynamicgut.run_simulation(range(3), single_models, pair_models, diet_file, density_file, test_folder)

    def test_bad_density_path(self, data_folder, model_files, pair_model_files, test_folder):
        single_models = [join(data_folder, m) for m in model_files]
        pair_models = [join(data_folder, m) for m in pair_model_files]
        diet_file = join(data_folder, 'initial_unhealthy_diet.json')
        density_file = join(data_folder, 'bad.csv')
        with pytest.raises(IOError):
            dynamicgut.run_simulation(range(3), single_models, pair_models, diet_file, density_file, test_folder)

    def test_bad_folder_path(self, data_folder, model_files, pair_model_files):
        single_models = [join(data_folder, m) for m in model_files]
        pair_models = [join(data_folder, m) for m in pair_model_files]
        diet_file = join(data_folder, 'initial_unhealthy_diet.json')
        density_file = join(data_folder, 'initial_density.csv')
        with pytest.raises(OSError):
            dynamicgut.run_simulation(range(3), single_models, pair_models, diet_file, density_file, '/bad_path')
