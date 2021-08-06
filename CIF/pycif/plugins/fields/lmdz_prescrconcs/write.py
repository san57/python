import xarray

from pycif.utils.classes.fluxes import Fluxes


def write(self, name, ic_file, flx, mode="a"):
    """Write prescribed species files for LMDZ

    Args:
        self (Fluxes): the Fluxes plugin
        ic_file (str): the file where to write fluxes
        flx (xarray.DataArray): fluxes data to write
        mode (str): 'w' to overwrite, 'a' to append
        """

    print("NO YET CODED")
