import pytest


def test_all():
    """ Run all of the tests. """
    return pytest.main(['--pyargs', 'dynamicgut', '-v']) == 0
