
from .conftests.adjtltest import dummy_config_adjtltest


def test_integration_adjtltest(dummy_config_adjtltest):
    """
    Integration test that runs the dummy_forward model.
    """

    # Test of the adjoint passes only if eps < 100
    assert dummy_config_adjtltest < 100
    
