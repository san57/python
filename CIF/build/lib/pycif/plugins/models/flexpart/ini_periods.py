import pandas as pd
import itertools
from pycif.utils.dates import date_range


def ini_periods(self, **kwargs):

    # Hardcoded 1h for now 
    datei = self.datei
    datef = self.datef
    period = getattr(self, "period", "1MS")
    subperiod = getattr(self, "subperiod", "")

    # List of sub-simulation windows
    self.subsimu_dates = date_range(
        datei,
        datef,
        period=period)

    # Time steps per day
    freq = '{:.0f}s'.format(3600)

    # List of time steps
    self.tstep_dates = \
                       {ddi: date_range(ddi, ddf, period=freq)
                        for ddi, ddf in zip(self.subsimu_dates[:-1],
                                            self.subsimu_dates[1:])}

    # List of input dates
    self.input_dates = \
                       {ddi: date_range(ddi, ddf,
                                        period=period,
                                        subperiod=subperiod)
                        for ddi, ddf in zip(self.subsimu_dates[:-1],
                                            self.subsimu_dates[1:])}

    # Merged list of time steps
    self.tstep_all = date_range(datei, datef,
                                period=self.periods, subperiod=freq)

