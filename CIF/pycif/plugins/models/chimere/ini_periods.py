from __future__ import division

import datetime

import numpy as np
import pandas as pd
from builtins import range

from pycif.utils.dates import date_range
from pycif.utils.netcdf import readnc


def ini_periods(self, **kwargs):
    datei = self.datei
    datef = self.datef
    nho = self.nho
    nhours = self.nhours

    # List of sub-simulation windows
    self.subsimu_dates = date_range(datei, datef, period=self.periods)

    # Time steps defined by the user
    nphour_ref = self.nphour_ref

    # Numbers of time steps per hour really used in the simulation
    self.tstep_dates = {}
    self.tstep_all = []
    self.nhour = []
    self.subtstep = []
    self.input_dates = {}

    for dd in self.subsimu_dates[:-1]:
        # time-steps in METEO.nc, computed by diagmet
        sdc = dd.strftime("%Y%m%d%H")
        met = "{}/METEO.{}.{}.nc".format(self.meteo.dir, sdc, nho)
        nbstep = readnc(met, ["nphourm"]).astype(int)

        # Loop on hours and check CFL
        self.tstep_dates[dd] = []
        for nh in range(nhours):
            nphour = nbstep[nh] if nphour_ref < nbstep[nh] else nphour_ref

            # Saving substep indexes for matching with observation
            self.subtstep.extend(list(range(1, nphour + 1)))
            self.nhour.extend(nphour * [nh + 1])

            # Frequency in seconds
            # TODO: Check with FORTRAN: nphour rounding?
            freq = "{}s".format(int(3600 // nphour))

            # List of time steps
            # TODO: what about chemical time steps?
            # the time step to really use is nphour*ichemstep
            ddhi = dd + datetime.timedelta(hours=nh)
            ddhe = dd + datetime.timedelta(hours=nh + 1)
            drange = list(pd.date_range(ddhi, ddhe, freq=freq).to_pydatetime())
            self.tstep_dates[dd].extend(drange[:-1])
            self.tstep_all.extend(drange[:-1])

        # List of dates for which inputs are needed
        self.input_dates[dd] = pd.date_range(
            dd, periods=nhours + 1, freq="1H"
        ).to_pydatetime()

        # Include last time step
        self.tstep_dates[dd].append(ddhe)
        self.tstep_dates[dd] = np.array(self.tstep_dates[dd])

    # Include very last time step
    self.tstep_all.append(ddhe)
    self.tstep_all = np.array(self.tstep_all)
