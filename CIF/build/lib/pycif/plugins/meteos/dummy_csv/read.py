#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pandas as pd


def read(meteo, name, tracdir, tracfile, varnames, dates, **kwargs):
    """Reads meteorology and links to the working directory

    Args:
        meteo (dictionary): dictionary defining the domain. Should include
        dirmeteo to be able to read the meteorology
        datei (datetime.datetime): initial date for the inversion window
        datef (datetime.datetime): end date for the inversion window
        workdir (str): path to the working directory where meteo files
                       should be copied
        logfile (str): path to the log file
        **kwargs (dictionary): extra arguments

    Return:
        None

    """

    # Read the reference meteo file
    ds = pd.read_csv(tracfile, index_col="date", parse_dates=True)

    # Putting data in the meteo Plugin
    meteo.data = ds

    return ds
