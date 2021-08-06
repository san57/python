from .conftests.variat import dummy_config_variat


def test_integration_variat(dummy_config_variat):
    """
    Integration test that runs the dummy_forward model.
    """
    
    print(dummy_config_variat)
    # Run the dummy model
    # result = Setup.run_simu({"def_file": dummy_config_file})
    #
    # Test of the adjoint passes only if eps < 100
    # assert result < 100
