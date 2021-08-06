import os

import numpy as np
import pandas as pd
import xarray

from pycif.utils.classes.fluxes import Fluxes
from pycif.utils.netcdf import save_nc


def write(self, name, ic_file, flx, mode="a"):
    """Write flux to AEMISSION CHIMERE compatible files.

    Args:
        self (Fluxes): the Fluxes plugin
        ic_file (str): the file where to write fluxes
        flx (xarray.DataArray): fluxes data to write
        mode (str): 'w' to overwrite, 'a' to append
        """

    # If mode is 'a' but file does not exit, switch to mode 'w'
    if mode == "a" and not os.path.isfile(ic_file):
        mode = "w"

    # Array shape
    nhours, nlev, nmerid, nzonal = np.shape(flx)

    # Dimensions
    spstrlen = 23
    dimnames = [
        "Time",
        "south_north",
        "west_east",
        "bottom_top",
        "SpStrLen",
        "DateStrLen",
        "Species",
    ]
    dimlens = [None, nmerid, nzonal, nlev, spstrlen, 19, 1]

    # Variables names, dimension and attributes
    varnames = ["Times", "species", "lon", "lat", name]
    vardims = [
        ("Time", "DateStrLen"),
        ("Species", "SpStrLen"),
        ("south_north", "west_east"),
        ("south_north", "west_east"),
        ("Time", "bottom_top", "south_north", "west_east"),
    ]
    dtypes = ["c", "c", "f", "f", "d"]
    units = ["", "", "degrees_east", "degrees_north", "molecule/cm2/s"]
    attributes = [
        {},
        {},
        {"long_name": "Longitude"},
        {"long_name": "Latitude"},
        {"long_name": "{} emissions".format(name)},
    ]

    # Variables to save
    times = [
        list(pd.to_datetime(d).strftime("%Y-%m-%d_%H:00:00"))
        for d in flx.time.values
    ]
    specs = [list(name.ljust(spstrlen))]
    lon = self.domain.zlon
    lat = self.domain.zlat

    variables = [times, specs, lon, lat, flx]

    save_nc(
        ic_file,
        variables,
        varnames,
        vardims,
        dimnames,
        dimlens,
        units,
        dtypes,
        mode=mode,
        attributes=attributes,
        format="NETCDF3_CLASSIC",
    )
