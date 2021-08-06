import datetime

import numpy as np


def dateslice(tracer, di, df):
    """Gets the temporal chunk of the control vector corresponding to a
    simulation period.

    Special cases:
        - one period spanning the whole inversion window:
            returns the whole table
        - initial conditions

    Args:
        tracer: the Tracer class with all its attributes
        di, df: start/end dates of the simulations period

    Returns:
        list of date indexes and list of dates (including the end date)

    """
    ddi = min(di, df)
    ddf = max(di, df)

    if hasattr(tracer, "dates") and hasattr(tracer, "period"):
        dates = np.append(
            tracer.dates, tracer.dates[-1] + datetime.timedelta(1)
        )
        chunk = np.where((dates[1:] > ddi) & (dates[:-1] <= ddf))[0]

    # When no period is specified and/or no date,
    # just take all the first element
    else:
        chunk = [0]

    return chunk
