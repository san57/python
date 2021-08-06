import numpy as np
import xarray as xr
from builtins import zip


def read(
    self,
    name,
    tracdir,
    tracfile,
    varname,
    dates,
    interpol_flx=False,
    tracer=None,
    **kwargs
):
    """Get fluxes from pre-computed fluxes and load them into a pycif
    variables

    Args:
        self: the model Plugin
        name: the name of the component
        tracdir, tracfile: flux directory and file format
        dates: list of dates to extract
        interpol_flx (bool): if True, interpolates fluxes at time t from
        values of surrounding available files

    """

    list_file_flx = [dd.strftime(tracfile) for dd in dates]

    # Reading fluxes for periods within the simulation window
    trcr_flx = []
    for dd, file_flx in zip(dates, list_file_flx):
        nc = xr.open_dataset(
            "{}/{}".format(tracdir, file_flx), decode_times=False
        )

        nlon = self.domain.nlon
        nlat = self.domain.nlat

        # Vector to map
        # Deals with polar boxes by sub-dividing them zonally
        # Also loops zonally for consistency with other call to gridded values
        flx = nc[varname].values
        flx0 = flx[:, 0]
        flx1 = flx[:, -1]

        flx = flx[:, 1:-1].reshape((-1, nlat - 2, nlon - 1))

        flx = np.append(
            flx,
            flx1[:, np.newaxis, np.newaxis] * np.ones((1, 1, nlon - 1)),
            axis=1,
        )
        flx = np.append(
            flx0[:, np.newaxis, np.newaxis] * np.ones((1, 1, nlon - 1)),
            flx,
            axis=1,
        )
        flx = np.append(flx, flx[:, :, np.newaxis, 0], axis=2)

        # Keeps only values for the corresponding month
        # Assumes monthly resolution
        if nc.dims["time"] == 12:
            month = dd.month
            flx = flx[month - 1]
        else:
            flx = flx[0]

        trcr_flx.append(flx)

    # Interpolating fluxes temporally between file values
    if interpol_flx:
        weights = []
        weights_inds = []
        for flx_file, flx, dd in zip(list_file_flx, trcr_flx, dates):
            inds = [
                k for k, flxx in enumerate(trcr_flx) if np.all(flx == flxx)
            ]
            w0 = dd - dates[inds[0]]
            w1 = dates[min(inds[-1] + 1, len(dates) - 1)] - dd
            dt = w1 + w0
            w0 = w0.total_seconds() / float(dt.total_seconds())
            w1 = w1.total_seconds() / float(dt.total_seconds())
            weights.append((w0, w1))
            weights_inds.append((inds[0], min(inds[-1] + 1, len(dates) - 1)))

        trcr_flx_interp = []
        for k, ((w0, w1), (i0, i1)) in enumerate(zip(weights, weights_inds)):
            trcr_flx_interp.append(trcr_flx[i0] * w1 + trcr_flx[i1] * w0)
        trcr_flx = trcr_flx_interp

    xmod = xr.DataArray(
        np.array(trcr_flx)[:, np.newaxis, ...],
        coords={"time": dates},
        dims=("time", "lev", "lat", "lon"),
    )

    return xmod
