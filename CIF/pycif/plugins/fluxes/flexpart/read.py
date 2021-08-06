import os
import pandas as pd
import pycif.utils.check as check
from pycif.utils import path
import xarray as xr
import datetime
import numpy as np
from pycif.utils.netcdf import readnc
from pycif.utils.dataarrays.reindex import reindex

# from pycif.utils.dates import j2d


def read(self, name, tracdir, tracfile, varnames, dates,
         interpol_flx=False, **kwargs):
    """Get fluxes from pre-computed fluxes and load them into a pycif
    variables

    Args:
        self: the model Plugin
        name: the name of the component
        tracdir, tracfic: flux directory and file format
        dates: list of dates to extract
        interpol_flx (bool): if True, interpolates fluxes at time t from
        values of surrounding available files

    """

    # Todo: can probably just get years in range dates[0], dates[-1]
    list_fic_flx = np.unique([dd.strftime(tracfile) for dd in dates])

    # Reading fluxes for periods within the simulation window
    trcr_flx = np.empty((0, *self.domain.zlat.shape), dtype=np.float)
    times = []
    for file_flx in list(list_fic_flx):
        # First load inside domain
        data, lat, lon, time_jd = \
            readnc(os.path.join(tracdir, file_flx),
                   [self.varname_flx, self.latname_flx,
                    self.lonname_flx, self.timename_flx])

        # Convert julian day (since 1-1-1900) to datetime
        times.extend(
            [datetime.datetime(1900, 1, 1) + datetime.timedelta(int(t))
             for t in time_jd]
        )

        # Convert to ng/m2/s
        numscale = np.float(getattr(self, 'numscale', 1.E12))
        data *= numscale / 3600.

        # Extract only data covering the inversion region
        ix0 = np.argmin(np.abs(lon - self.domain.lon_in[0]))
        iy0 = np.argmin(np.abs(lat - self.domain.lat_in[0]))

        flx_reg_in = data[:, iy0:iy0 + self.domain.nlat,
                     ix0:ix0 + self.domain.nlon]
        flx_reg_in = flx_reg_in.reshape(flx_reg_in.shape[0], -1)
        
        # Loading outside data
        out_file = os.path.join(
            os.path.dirname(os.path.join(tracdir, file_flx)),
            times[0].strftime(self.file_glob)
        )
        data, lat, lon, time_jd = \
            readnc(out_file,
                   [self.varname_flx, self.latname_flx,
                    self.lonname_flx, self.timename_flx])

        # Extract data outside nest domain
        flx_reg_out = \
            np.delete(data.reshape(data.shape[0], -1),
                      self.domain.raveled_indexes_glob, 1)
        
        # Concatenate
        trcr_flx = np.append(
            trcr_flx,
            np.append(flx_reg_in, flx_reg_out, axis=1)[:, :, np.newaxis],
            axis=0)
    
    # Put data into dataarry
    xmod = xr.DataArray(trcr_flx[:, np.newaxis],
                        coords={'time': np.array(times)},
                        dims=('time', 'lev', 'lat', 'lon'))
    
    # Reindex to required dates
    xmod = reindex(xmod, levels={"time": dates.astype(np.datetime64)})

    #
    #
    # # TODO: take care if several files are read
    # # TODO: scale flux contribution by area weight for boxes
    # # TODO: consider storing fluxes at original time resolution and
    # #       interpolate as needed
    #
    # flx = np.ndarray((self.ndates, self.domain.nlat, self.domain.nlon))
    #
    # # Interpolate fluxes to start time of control period
    # for ddt in range(self.ndates):
    #     if interpol_flx:
    #         flx[ddt, :, :] = xmod.interp(time=self.dates[ddt])
    #     else:
    #         flx[ddt, :, :] = xmod.sel(time=self.dates[ddt], method='nearest')

    return xmod
