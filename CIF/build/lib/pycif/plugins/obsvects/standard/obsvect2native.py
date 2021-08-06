from __future__ import absolute_import

import pandas as pd


def obsvect2native(self, datei, datef, mode, runsubdir, workdir):
    """Convert elements of the observation vector to the native model
    resolution.
    At the moment, it only includes arithmetic temporal means of localized
    observations.

    In the future, super observations will be processed here.

    Args:
        self (pycif.utils.classes.obsvects.ObsVect): observation vector object
        datei (datetime.datetime): Start date of current sub-simulations
        datef (datetime.datetime): End date
        mode (str): running mode of the model; one of 'fwd', 'tl', 'adj'
        runsubdir (str): sub-simulation directory
        workdir (str): pycif work directory

    Returns:
        pandas.DataFrame
    """

    # Reversing dates if adjoint run
    ddi = min(datei, datef)
    ddf = max(datei, datef)

    datastore = self.datastore

    # If empty datastore, do nothing
    if datastore.size == 0:
        return datastore

    # Get observations for the sub-period
    mask = (
        datastore.index + pd.to_timedelta(datastore["duration"], unit="h")
        > ddi
    ) & (datastore.index < ddf)
    y0 = datastore.loc[mask]

    return y0
