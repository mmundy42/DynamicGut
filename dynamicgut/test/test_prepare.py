from os.path import join, basename
from shutil import rmtree
from micom import load_pickle
import pytest

import dynamicgut


class TestPrepare:

    def test_serial(self, data_folder, model_files, test_folder):
        pair_folder = join(test_folder, 'pair_models')
        single_models = [join(data_folder, m) for m in model_files]
        pair_models = dynamicgut.prepare(single_models, pair_folder)
        assert len(pair_models) == 3
        for pair in pair_models:
            test_model = load_pickle(pair)
            good_model = load_pickle(join(data_folder, basename(pair)))
            assert len(test_model.reactions) == len(good_model.reactions)
            assert len(test_model.metabolites) == len(good_model.metabolites)
        rmtree(pair_folder)

    def test_multi(self, data_folder, model_files, test_folder):
        pair_folder = join(test_folder, 'pair_models')
        single_models = [join(data_folder, m) for m in model_files]
        pair_models = dynamicgut.prepare(single_models, pair_folder, n_processes=2)
        assert len(pair_models) == 3
        for pair in pair_models:
            test_model = load_pickle(pair)
            good_model = load_pickle(join(data_folder, basename(pair)))
            assert len(test_model.reactions) == len(good_model.reactions)
            assert len(test_model.metabolites) == len(good_model.metabolites)
        rmtree(pair_folder)

    def test_optimize(self, data_folder, model_files, test_folder):
        pair_folder = join(test_folder, 'pair_models')
        single_models = [join(data_folder, m) for m in model_files]
        pair_models = dynamicgut.prepare(single_models, pair_folder, optimize=True, n_processes=2)
        assert len(pair_models) == 3
        rmtree(pair_folder)

    def test_solver(self, data_folder, model_files, test_folder):
        pair_folder = join(test_folder, 'pair_models')
        single_models = [join(data_folder, m) for m in model_files]
        pair_models = dynamicgut.prepare(single_models, pair_folder, solver='glpk')
        assert len(pair_models) == 3
        rmtree(pair_folder)

    def test_bad_single_path(self, data_folder, model_files, test_folder):
        pair_folder = join(test_folder, 'pair_models')
        single_models = [join(data_folder, m) for m in model_files]
        single_models.append(join(data_folder, 'bad.xml'))
        with pytest.raises(IOError):
            dynamicgut.prepare(single_models, pair_folder)

    def test_bad_folder_path(self, data_folder, model_files):
        single_models = [join(data_folder, m) for m in model_files]
        with pytest.raises(OSError):
            dynamicgut.prepare(single_models, '/bad_path')
