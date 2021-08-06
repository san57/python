import datetime

import numpy as np
import pandas as pd
from builtins import zip


def date2datetime(date, time=(0, 0)):
    """Converts a date to a datetime assuming that if no hour is specified,
    00:00 is the standard

    Args:
        date (date): date to convert
        time (tuple): time as a tuple of integers (hour, minute)

    Return:
        datetime
    """

    if type(date) is datetime.datetime:
        return date

    else:
        return datetime.datetime.combine(date, datetime.time(time[0], time[1]))


def date_range(datei, datef, period="1000000H", subperiod=""):
    """Generate a range of dates between two dates, with the specified time
    interval in days. Fills gaps in the end of months to make sure that the
    first day of each month is returned and that a sub-period is never
    shorter than $period

    Args:
        datei (datetime): starting date
        datef (datetime): end date
        period (str): pandas frequency
        subperiod (str): pandas sub-frequency

    """

    # If the period is empty, keeps the full period as a sub-period
    if period == "":
        list_dates = [datei]
        return np.array(list_dates)

    # Otherwise, split simulation window along base periods
    drange = pd.date_range(datei, datef, freq=period)

    # If period is longer than the length of the simulation window,
    # returns a very far away end date
    # Otherwise, just put the end date at the end of the list
    # This is done in case the last time stamp should be use
    # for interpolation in the last hours of the run
    if datei + pd.to_timedelta(period) > datef:
        return (
            drange.union([drange[-1] + pd.to_timedelta(period)])
            .to_series()
            .dt.to_pydatetime()
        )
    else:
        drange = drange.union([datef])

    # One can specify sub-periods to split main periods
    if subperiod == "":
        list_dates = drange.to_series().dt.to_pydatetime()

    else:
        list_dates = []
        for d0, d1 in zip(drange[:-1], drange[1:]):
            drange = pd.date_range(d0, d1, freq=subperiod.format(period))
            list_dates.extend(drange.to_series().dt.to_pydatetime())

        if datef not in list_dates:
            list_dates.append(datef)

        for d0, d1 in zip(np.array(list_dates)[:-1], np.array(list_dates)[1:]):
            if d1 - d0 < pd.to_timedelta(drange.freq) / 2.0:
                list_dates.remove(d0)

    return np.array(list_dates)
