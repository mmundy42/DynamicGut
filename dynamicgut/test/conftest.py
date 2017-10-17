import pytest
from os.path import join, abspath, dirname
from tempfile import mkdtemp


@pytest.fixture(scope='session')
def data_folder():
    dynamicgut_folder = abspath(join(dirname(abspath(__file__)), '..'))
    return join(dynamicgut_folder, 'test', 'data')


@pytest.fixture(scope='session')
def test_folder():
    return mkdtemp(prefix='dynamicgut')


@pytest.fixture(scope='session')
def model_files():
    return ['Bacteroides_thetaiotaomicron_VPI_5482.xml',
            'Escherichia_coli_str_K_12_substr_MG1655.xml',
            'Eubacterium_rectale_ATCC_33656.xml']
