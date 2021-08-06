import numpy as np


def write(self, name, flx_file, flx, mode="w"):
    """Write fluxes for the dummy_txt Gaussian model.

    Args:
        self (Fluxes): the Fluxes plugin
        file (str): the file where to write fluxes
        flx_fwd, flx_tl: fluxes data to write
        """

    np.savetxt(flx_file, flx, delimiter=",")
