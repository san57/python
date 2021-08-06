import os
import pandas as pd
import pycif.utils.check as check
from pycif.utils import path
import xarray as xr
import datetime
import numpy as np
from pycif.utils.netcdf import readnc
# from pycif.utils.dates import j2d


def read_glob(self, name, tracdir, tracfic, dates,
                interpol_flx=False, **kwargs):
    """Get global fluxes from pre-computed fluxes and load them into a pycif
    variables

    Args:
        self: the model Plugin
        name: the name of the component
        tracdir, tracfic: flux directory and file format
        dates: list of dates to extract
        interpol_flx (bool): if True, interpolates fluxes at time t from
        values of surrounding available files

    Note:
        This was originally copied from ../flexpart/read.py. May eventually 
        be moved to a different plugin

    """

    # Available files in the directory
    list_files = os.listdir(tracdir)
    dates_available = []
    for fic in list_files:
        try:
            dates_available.append(datetime.datetime.strptime(fic, tracfic))
        except:
            continue

    dates_available = np.array(dates_available)
    
    # Todo: can probably just get years in range dates[0], dates[-1]
    list_fic_flx = np.unique([dd.strftime(tracfic) for dd in dates])
    
    # Reading fluxes for periods within the simulation window
    trcr_flx = []
    times = []
    #    for dd, fic_flx in zip(dates, list_fic_flx):
    trcr_flx = []
    for file_flx in list(list_fic_flx):
        data, lat, lon, time_jd =  readnc(os.path.join(tracdir, file_flx),
                                       [self.varname_flx, self.latname_flx,
                                        self.lonname_flx, self.timename_flx])

        # Convert julian day (since 1-1-1900) to datetime
        for t in time_jd:
            times.append(datetime.datetime(1900,1,1) + datetime.timedelta(int(t)))

        # Convert to ng/m2/s
        numscale = np.float(getattr(self, 'numscale', 1.E12) )
        data *= numscale/3600.
        
        trcr_flx.append(data[:, :, :])
        
    xmod = xr.DataArray(trcr_flx[0],
                        coords={'time': times},
                        dims=('time', 'lat', 'lon'))



    # TODO: take care if several files are read
    # TODO: scale flux contribution by area weight for boxes
    # TODO: consider storing fluxes at original time resolution and
    #       interpolate as needed

    flx = np.ndarray((self.ndates, self.domain.nlat_glob, self.domain.nlon_glob))

    # Interpolate fluxes to start time of control period
    for ddt in range(self.ndates):
        if interpol_flx:
            flx[ddt, :, :] = xmod.interp(time=self.dates[ddt])
        else:
            flx[ddt, :, :] = xmod.sel(time=self.dates[ddt], method='nearest')
    
    return flx
