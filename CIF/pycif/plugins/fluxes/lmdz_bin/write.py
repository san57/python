import numpy as np
import xarray

from pycif.utils.classes.fluxes import Fluxes


def write(self, name, flx_file, flx, mode="a"):
    """Write flux to LMDZ-DISPERSION compatible files.
    The shape follows the LMDZ physical vectorial shape grid.
    Each line of the output binary file includes.

    Args:
        self (Fluxes): the Fluxes plugin
        flx_file (str): the file where to write fluxes
        flx (xarray.Dataset): fluxes data to write
        """

    flx_fwd = flx["fwd"].values[:, 0]
    flx_tl = flx["tl"].values[:, 0]
    np.transpose([flx_fwd, flx_tl], axes=(0, 3, 2, 1)).T.tofile(flx_file)
