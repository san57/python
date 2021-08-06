import datetime

import pandas as pd

from pycif.utils.check import info


def crop_monitor(datastore, datei, datef, **kwargs):
    """Crops observation datasets to keep observations whose duration fits
    entirely during the simulation period

    Args:
        datastore (pd.DataFrame): observation dataset
        datei (datetime.datetime): start date
        datef (datetime.datetime): end date

    Returns:
        pd.DataFrame: Cropped dataframe
    """

    info("Cropping obsvect.datastore to simulation window")
    info("{} to {}".format(datei, datef))

    mask = (datastore.index >= datei) & (
        datastore.index + pd.to_timedelta(datastore["duration"], unit="h")
        <= datef
    )
    return datastore.loc[mask]
