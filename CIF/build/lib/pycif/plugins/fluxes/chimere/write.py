import os

from netCDF4 import Dataset
import numpy as np
import pandas as pd
import xarray

from pycif.utils.classes.fluxes import Fluxes
from pycif.utils.netcdf import save_nc


def write(self, name, flx_file, flx, mode="a"):
    """Write flux to AEMISSION CHIMERE compatible files.

    Args:
        self (Fluxes): the Fluxes plugin
        flx_file (str): the file where to write fluxes
        flx (xarray.DataArray): fluxes data to write
        mode (str): 'w' to overwrite, 'a' to append
        """

    # If mode is 'a' but file does not exit, switch to mode 'w'
    if mode == "a" and not os.path.isfile(flx_file):
        mode = "w"

    # Variables from self
    lon = self.domain.zlon
    lat = self.domain.zlat

    nlev = self.nlevemis

    write_AEMISSIONS(name, flx_file, flx, mode, lon, lat, nlev)


def write_AEMISSIONS(
    name, flx_file, flx, mode="a", lon=None, lat=None, nlev=1
):
    """Auxiliary function that can be used outside pycif"""

    # Array shape
    nmerid, nzonal = np.shape(lon)

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
    # varnames = ['Times', 'species', 'lon', 'lat', name]
    varnames = [name]
    vardims = [("Time", "bottom_top", "south_north", "west_east")]
    dtypes = ["d"]
    units = ["molecule/cm2/s"]
    attributes = [{"long_name": "{} emissions".format(name)}]
    variables = [flx]

    # Add auxiliary units if not present
    times = [
        list(pd.to_datetime(d).strftime("%Y-%m-%d_%H:00:00"))
        for d in flx.time.values
    ]
    
    auxnames = ["Times", "lon", "lat"]
    auxdims = [
        ("Time", "DateStrLen"),
        ("south_north", "west_east"),
        ("south_north", "west_east"),
    ]
    auxtypes = ["c", "f", "f"]
    auxunits = ["", "degrees_east", "degrees_north"]
    auxattributes = [
        {},
        {"long_name": "Longitude"},
        {"long_name": "Latitude"}
    ]
    auxvariables = [times, lon, lat]
    
    for vn, vd, dt, un, at, var \
            in zip(auxnames, auxdims, auxtypes, auxunits, auxattributes,
                   auxvariables):
        if mode == "w":
            continue

        with Dataset(flx_file, 'a') as f:
            if vn in f.variables:
                continue
            
            varnames.append(vn)
            vardims.append(vd)
            dtypes.append(dt)
            units.append(un)
            attributes.append(at)
            variables.append(var)
    
    # Variables to save

    # Output data
    ljust_specs_in = [list(name.ljust(spstrlen))]
    specs_in = [name]

    # Append species to existing file
    # add_species = False
    # if os.path.isfile(flx_file) and mode == 'a':
    #     with Dataset(flx_file, 'a') as f:
    #         ljust_specs_in = f.variables['species'][:]
    #         specs_in = [''.join(p).strip() for p in ljust_specs_in]
    #
    #     if name not in specs_in:
    #         add_species = True
    #         specs_in = specs_in + [name]
    #         ljust_specs_in = [list(s.ljust(spstrlen)) for s in specs_in]

    # if add_species:
    #     with Dataset(flx_file, 'a') as f:
    #         print(__file__)
    #         import code
    #         code.interact(local=dict(locals(), **globals()))
    #

    # variables = [times, ljust_specs_in, lon, lat, flx]
    try:
        save_nc(
            flx_file,
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
    except:
        print(__file__)
        import code
        code.interact(local=dict(locals(), **globals()))