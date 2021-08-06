import numpy as np
import xarray as xr

from .utils import grib_file_reader


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
    # tracfile can be a list of same lenght as dates
    try:
        if len(tracfile) != len(dates):
            raise Exception(
                "Try read EDGAR files from a list of dates and a "
                "list of files, but not of same length:\n{}\n{}".format(
                    tracfile, dates
                )
            )
        list_files = tracfile[:]
    except TypeError:
        list_files = len(dates) * [tracfile]

    xout = []
    for dd, dd_file in zip(dates, list_files):
        dir_dd = dd.strftime(tracdir)
        spec = grib_file_reader("{}/{}".format(dir_dd, dd_file), varnames)
        xout.append(spec)

    # Flipping upside down to follow convention from other models
    xmod = xr.DataArray(
        np.array(xout)[:, ::-1],
        coords={"time": dates},
        dims=("time", "lev", "lat", "lon"),
    )

    return xmod
