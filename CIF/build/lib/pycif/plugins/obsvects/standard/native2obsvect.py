from __future__ import absolute_import

import pandas as pd


def native2obsvect(self, datastore, datei, datef, runsubdir, workdir):
    """Converts simulations from the native resolution of the model
    to the resolution of the control vector.
    Currently it includes temporal arithmetic averaging only.
    In the future, super observations should be processed here

    Args:
        self (pycif.utils.classes.obsvects.ObsVect): the observation vector
        datastore (pandas.DataFrame): the simulations extracted from the model
        datei (datetime.datetime): starting date for the sub-simulation to
        process
        datef (datetime.datetime): end date
        runsubdir (str): path to the current simulation directory
        workdir (str): path to the pycif working directory

    """

    di = min(datei, datef)
    df = max(datei, datef)

    mask = (
        self.datastore.index
        + pd.to_timedelta(self.datastore["duration"], unit="h")
        > di
    ) & (self.datastore.index < df)

    # Re-aggregating into self.datastore
    self.datastore.loc[mask, ["sim", "sim_tl"]] += datastore.loc[
        :, ["sim", "sim_tl"]
    ].values
