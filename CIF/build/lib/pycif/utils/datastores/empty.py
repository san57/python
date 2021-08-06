import pandas as pd


def init_empty():
    """Generates a empty pyCIF datastores with correct columns"""

    columns = [
        "station",
        "network",
        "parameter",
        "lon",
        "lat",
        "alt",
        "i",
        "j",
        "level",
        "obs",
        "obserror",
        "sim",
        "sim_tl",
        "tstep",
        "tstep_glo",
        "dtstep",
        "duration",
    ]

    return pd.DataFrame({c: [] for c in columns}, index=pd.to_datetime([]))
