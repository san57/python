from __future__ import division

import datetime
import glob
import os

import numpy as np
import xarray as xr

from pycif.utils.check import info


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

    # tracfile can be a list of same length as dates
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

    # Reading fluxes for periods within the simulation window
    trcr_flx = []
    for dd, dd_file in zip(dates, list_files):
        file_flx = dd.strftime(dd_file)
        dir_flx = dd.strftime(tracdir)

        if not os.path.isfile("{}/{}".format(dir_flx, file_flx)) and getattr(
            self, "closest_year", False
        ):
            info(
                "Warning: could not find correct year for EDGAR; "
                "using closest available one"
            )
            list_dates = [
                datetime.datetime.strptime(os.path.basename(f), tracfile)
                for f in glob.glob("{}/v50_*nc".format(dir_flx))
            ]
            delta_dates = np.abs(dd - np.array(list_dates))
            file_flx = list_dates[np.argmin(delta_dates)].strftime(tracfile)

        nc = xr.open_dataset(
            "{}/{}".format(dir_flx, file_flx), decode_times=False
        )
        trcr_flx.append(nc[varnames].values)

    xmod = xr.DataArray(
        np.array(trcr_flx)[:, np.newaxis, ...],
        coords={"time": dates},
        dims=("time", "lev", "lat", "lon"),
    )

    return xmod
