import datetime
import os

import numpy as np
import xarray as xr
from netCDF4 import Dataset

from pycif.utils.netcdf import readnc


def read(
    self,
    name,
    tracdir,
    tracfile,
    varnames,
    dates,
    interpol_flx=False,
    comp_type=None,
    **kwargs
):
    """Get fluxes from pre-computed fluxes and load them into a pyCIF
    variables

    Args:
        self: the fluxes Plugin
        name: the name of the component
        tracdir, tracfile: flux directory and file format
        dates: list of dates to extract
        interpol_flx (bool): if True, interpolates fluxes at time t from
        values of surrounding available files

    """
    # Check the type of limit condition to check
    if comp_type is None:
        raise Exception(
            "Trying to read limit conditions for CHIMERE, "
            "but did not specify the type"
        )

    # Read INI_CONCS
    if comp_type == "inicond":
        ic_file = min(dates).strftime("{}/{}".format(tracdir, tracfile))
        with Dataset(ic_file, "r") as f:
            data = f.variables[name][:][np.newaxis, :]
            xmod = xr.DataArray(
                data,
                coords={"time": [min(dates)]},
                dims=("time", "lev", "lat", "lon"),
            )

    # Read Lateral boundary conditions
    elif comp_type in ["latcond", "topcond"]:
        # Available files in the directory
        list_files = os.listdir(tracdir)
        list_available = []
        for bc_file in list_files:
            try:
                list_available.append(
                    datetime.datetime.strptime(bc_file, tracfile)
                )
            except BaseException:
                continue

        list_available = np.array(list_available)
        list_available.sort()

        # Reading required fluxes files
        trcr_bc = []
        for dd in dates:
            delta = dd - list_available
            mask = delta >= datetime.timedelta(0)
            imin = np.argmin(delta[mask])
            fdates = list_available[mask][imin]

            # Getting the data
            filein = fdates.strftime("{}/{}".format(tracdir, tracfile))
            spec = "lat_conc" if comp_type == "latcond" else "top_conc"
            data, times, specs = readnc(filein, [spec, "Times", "species"])

            # Get the correct date and species index
            ispec = ["".join(c).strip() for c in specs].index(name)
            idate = [
                datetime.datetime.strptime("".join(d), "%Y-%m-%d_%H:%M:%S")
                for d in times
            ].index(dd)

            # Appending
            trcr_bc.append(data[idate, ..., ispec])

        # Putting the data into an xarray
        # Adding an empty latitude axis
        if comp_type == "latcond":
            xout = np.array(trcr_bc)[..., np.newaxis, :]
        else:
            xout = np.array(trcr_bc)[:, np.newaxis, ...]

        xmod = xr.DataArray(
            xout, coords={"time": dates}, dims=("time", "lev", "lat", "lon")
        )
    else:
        raise Exception(
            "Could not recognize the type of boundary condition "
            "to read in CHIMERE: {}".format(comp_type)
        )

    return xmod
