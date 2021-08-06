import datetime
import os

import numpy as np
import xarray as xr

from pycif.utils.netcdf import readnc


def read(
    self,
    name,
    tracdir,
    tracfile,
    varnames,
    dates,
    interpol_flx=False,
    tracer=None,
    model=None,
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

    # Replace tracfile by available information from model
    if tracfile == "":
        tracfile = model.fluxes.file

    # Available files in the directory
    list_files = os.listdir(tracdir)
    list_available = []
    for flx_file in list_files:
        try:
            list_available.append(
                datetime.datetime.strptime(flx_file, tracfile)
            )
        except BaseException:
            continue

    list_available = np.array(list_available)

    # Reading required fluxes files
    trcr_flx = []
    for dd in dates:
        delta = dd - list_available
        mask = delta >= datetime.timedelta(0)
        imin = np.argmin(delta[mask])
        fdates = list_available[mask][imin]

        filein = fdates.strftime("{}/{}".format(tracdir, tracfile))

        data, times = readnc(filein, [name, "Times"])

        # Get the correct hour in the file
        times = [
            datetime.datetime.strptime(
                str(b"".join(s), "utf-8"), "%Y-%m-%d_%H:%M:%S"
            )
            for s in times
        ]
        hour = int((dd - times[0]).total_seconds() // 3600)

        trcr_flx.append(data[hour, ...])
    
    # Building a xarray
    xmod = xr.DataArray(
        trcr_flx, coords={"time": dates}, dims=("time", "lev", "lat", "lon")
    )

    return xmod
