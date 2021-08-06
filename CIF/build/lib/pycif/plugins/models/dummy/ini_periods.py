from builtins import zip

from pycif.utils.dates import date_range


def ini_periods(self, **kwargs):
    datei = self.datei
    datef = self.datef

    # List of sub-simulation windows
    self.subsimu_dates = date_range(datei, datef, period=self.periods)

    # Fixed time step (in hours): list of time steps for each sub-simulation
    self.tstep_dates = {
        ddi: date_range(ddi, ddf, period=self.tstep)
        for ddi, ddf in zip(self.subsimu_dates[:-1], self.subsimu_dates[1:])
    }

    # Fixed time step (in hours): list of required input time steps
    self.input_dates = {
        ddi: date_range(ddi, ddf, period=self.tstep)
        for ddi, ddf in zip(self.subsimu_dates[:-1], self.subsimu_dates[1:])
    }

    # Merged list of time steps
    self.tstep_all = date_range(
        datei, datef, period=self.periods, subperiod=self.tstep
    )
