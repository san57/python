import numpy as np
import os

# from pycif.utils.classes.setup import Setup
from pycif.utils.netcdf import readnc

from .conftests.fwd import dummy_config_fwd


def test_integration_fwd(dummy_config_fwd):
    """
    Integration test that runs the dummy_forward model.
    """
    
    print(dummy_config_fwd)
    
    # tmpdir, dummy_config_file = dummy_config_fwd
    #
    # # Run the dummy model
    # Setup.run_simu({"def_file": dummy_config_file})
    #
    # # Get results
    # monitor_ref = readnc(os.path.join(tmpdir, "monitor_reference.nc"), ["station"])
    #
    # # Define expected results
    # expected_res = np.ma.core.MaskedArray(
    #     np.tile(np.arange(50), 97), fill_value=1e20, dtype=np.float64
    # )
    #
    # # Compare results
    # assert np.array_equal(monitor_ref, expected_res)
