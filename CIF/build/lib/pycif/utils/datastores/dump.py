#!/usr/bin/env python
# -*- coding: utf-8 -*-
import datetime
import os

import numpy as np
import pandas as pd
import xarray as xr

from pycif.utils.check import info


def dump_datastore(
    datastore,
    file_monit=None,
    workdir="",
    dump=True,
    spec2dump=None,
    col2dump=[],
    dump_default=True,
    index=None,
    mode="a",
    dump_type="nc",
    nc_attributes=None,
    **kwargs
):
    """Dumps the datastore to a dedicated NetCDF file. If no file path is
    specified, put it into pyCIF working directory as default.

    The format is inspired from ObsPack data format.
    It includes the following variables (some are optional):
    - station
    - network
    - species/parameter
    - time_components (integer time components)
    - time (NetCDF compatible seconds since 1970)
    - obs_period
    - tstep (time step in the model)
    - dtstep (number of model time step corresponding to the observation
    duration)
    - longitude
    - latitude
    - i, j = model grid cell corresponding to the measurement
    - altitude of inlet (a.s.l)
    - elevation of inlet (a.g.l)
    - observation value
    - obserror
    - prior value
    - posterior value

    Args:
        datastore (pd.DataFrame): dataframe to dump
        dump (bool) : dumps the data if True. Default is True
        file_monit (str) : file where to dump the data
        workdir (str): working directory
        col2dump (list[str]): columns to dump
        spec2dump (list[str]): species to dump
        index (str): either None to use the datastore index for dates,
                     or a column name
        mode (str): write or append to existing monitor
        dump_type (str): 'nc' for NetCDF or 'csv'
        nc_attributes (dict): global attributes to save in the NetCDF
        **kwargs (dictionary) : any additional options needed later

    """

    # Default columns to dump anyway
    if dump_default:
        col2dump = [
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
        ] + col2dump

    # If not directory is set, save to $workdir/monitor
    if file_monit == "" or file_monit is None:
        if workdir == "":
            info("Invalid output file. Dumping nothing")
            return
        else:
            file_monit = "{}/monitor.{}".format(workdir, dump_type)

    # If file exists and mode is 'a', reads existing monitor and merge
    if os.path.isfile(file_monit):
        if mode == "a":
            df_in = read_datastore(file_monit, dump_type=dump_type)
            datastore = pd.concat([datastore, df_in])
            datastore.drop_duplicates(inplace=True)

        os.remove(file_monit)

    if spec2dump is None:
        spec2dump = list(datastore.keys())

    # Dumping
    info("Dumping datastore to {}".format(file_monit))

    df0 = pd.DataFrame({}, index=datastore.index)
    for col in col2dump:
        if col in datastore.columns:
            df0[col] = datastore[col].values
        else:
            df0[col] = np.nan

    if dump_type == "csv":
        df0.to_csv(file_monit, index_label="date", date_format="%Y%m%d%H%M")

    elif dump_type == "nc":
        datastore.index.names = ["index"]
        ds = df0.to_xarray()

        # Saving attributes to the monitor
        history = datetime.datetime.now().strftime(
            "Created on %d-%m-%Y %H:%M:%S"
        )
        if nc_attributes is None:
            nc_attributes = {"history": history}

        else:
            nc_attributes["history"] = history

        for key in nc_attributes:
            ds.attrs[key] = nc_attributes[key]

        # Saving to NetCDF
        ds.to_netcdf(file_monit)

    else:
        raise TypeError(
            "pyCIF can dump observation files only as txt or "
            "nc files. {} was requested. Please check your "
            "set-up".format(dump_type)
        )


def read_datastore(
    file_monitor, chunksize=1e5, col2dump=[], dump_type="nc", **kwargs
):
    """Reads a datastore from a pre-defined text file

    The format must be:

    station_network_species year month day hour minute lon lat spec error alt

    Args:
        file_monitor (str): directory to the monitor file

    Kwargs:
        chunksize (int): size of chunks to read the monitor file. Avoids
            killing the memory with big observation files. Default is 1e5

    Returns:
        DataFrame with the observations

    """

    # Doesn't try to read if file_monitor is None
    if file_monitor is None:
        raise IOError

    if dump_type == "nc":
        ds = xr.open_dataset(file_monitor)

        # Converting strings
        for s in ds:
            if "|S" in str(ds[s].dtype):
                ds[s] = ds[s].astype(str)

        df = ds.to_dataframe()
        # df.index = df.index.to_frame()['index']
        setattr(df, "nc_attributes", ds.attrs)

    elif dump_type == "csv":
        # Reading the file by chunks in case it is a very big file
        df = pd.concat(
            [d for d in pd.read_csv(file_monitor, chunksize=chunksize)]
        )
        df["date"] = pd.to_datetime(df["date"], format="%Y%m%d%H%M")
        df.set_index("date", inplace=True)
        df.sort_index(inplace=True)

        if "sim_period" in df and not np.any(np.isnan(df["sim_period"])):
            df["sim_period"] = pd.to_datetime(
                df["sim_period"].astype(int), format="%Y%m%d%H%M"
            )

    # Re-ordering the dataframe before returning
    df.sort_index(inplace=True)

    return df
