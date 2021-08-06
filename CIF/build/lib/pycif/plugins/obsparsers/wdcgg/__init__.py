# -*- coding: utf-8 -*-
import os

import numpy as np
import pandas as pd
from builtins import map
from builtins import range
from builtins import zip

from pycif.utils.check import info
from .headers import get_header, parse_header
from .utils import parse_file_name, shiftdate, rescale, convert_unit


def do_parse(
    self,
    obs_file,
    maxlen=300,
    default_unit="ppm",
    default_tz="utc",
    default_duration=1,
    na_values=-999,
    err_na_values=-999,
    spec=None,
    extract=["obserror", "unit", "lon", "lat", "alt", "tz"],
    **kwargs
):
    """Parse function for a file from WDCGG

    Args:
        obs_file (str) :
            Path to input file
        maxlen (int):
            Maximum possible length for a WDCGG header. Default is `300`
        default_unit (str):
            Default unit to use for converting the processed species to ppm if
            not explicitly specified in the WDCGG file. Default is ppm
        default_tz (str):
            Default time zone to shift observations to UTC time if not
            explicitly specified in the WDCGG file. Default is utc
        default_duration (str):
            Default duration of measurements in hours. Default is 1.
        na_values (str):
            Pattern for NaN values. Default is -999
        spec (str, optional) :
            Species to extract. Default is None.
            If `MCF` is specified, it is translated to `CH3CCl3`.
        extract (list[str], optional):
            List of parameters to extract.
            Default is `[error, lon, lat, alt, tz]`
        logfile (str, optional) :
            File name for verbose entries

    Returns:
        pandas.DataFrame :
            Dataframe from input file df[parameter][station]

    """

    info("Reading observation file: {}".format(os.path.basename(obs_file)))

    # Get default unit from species description if available
    if hasattr(self, "default_unit"):
        default_unit = self.default_unit

    # Get default duration from species description if available
    if hasattr(self, "default_duration"):
        default_duration = self.default_duration

    # Get default na_values from species description if available
    if hasattr(self, "na_values"):
        na_values = self.na_values

    # Get default na_values from species description if available
    if hasattr(self, "err_na_values"):
        err_na_values = self.err_na_values

    # Scans file to get the header
    header = get_header(obs_file, maxlen)

    # Does not read empty files
    if len(header) == 0:
        info("{} is empty. Not reading it".format(obs_file))
        return pd.DataFrame({})

    else:
        # Get spec info either from the function argument
        # or from the file name
        file_infos = parse_file_name(obs_file)

        if spec is None:
            spec = file_infos["parameter"]

        list_extract = [spec] + extract

        # Get the content of columns and extra information from the header
        names, columns, date_ids, extra = parse_header(
            header, spec, list_extract, default_unit, default_tz
        )

        # Reads the file with Pandas
        df = pd.read_csv(
            obs_file,
            delim_whitespace=True,
            skiprows=len(header),
            usecols=date_ids + columns,
            parse_dates=[list(range(len(date_ids)))],
            infer_datetime_format=True,
            quoting=3,
            header=None,
            na_values=na_values,
        )

        # Rename columns according to standard names
        df.rename(columns=dict(list(zip(columns, names))), inplace=True)
        df.rename(columns={"_".join(map(str, date_ids)): "time"}, inplace=True)

        # Set the data frame index as time
        df.set_index("time", inplace=True)

        # Deals with problems in reading the dates
        # Some hour values at 99 and 24
        if not isinstance(df.index, pd.DatetimeIndex):
            index = df.index.values.astype("str")

            # Removes hours > 24 and minutes > 60
            hours = np.array(
                [ln.split(":")[0].split(" ")[-1] for ln in index]
            ).astype(int)
            df = df[hours <= 24]
            index = index[hours <= 24]

            minutes = np.array([ln.split(":")[1] for ln in index]).astype(int)
            df = df[minutes <= 60]
            index = index[minutes <= 60]

            # Shift hours at 24 by 1 hour back in time to be compatible for
            # reading
            # then forward
            mask24 = np.core.defchararray.find(index, "24:") > 0
            index = np.core.defchararray.replace(index, "24:", "23:")

            mask60 = np.core.defchararray.find(index, ":60") > 0
            index = np.core.defchararray.replace(index, ":60", ":00")

            index = pd.to_datetime(pd.Series(index))
            index[mask24] += pd.DateOffset(hours=1)
            index[mask60] += pd.DateOffset(hours=1)

            df.index = index

        # Shifting dates depending on time zone, then removing corresponding
        # key
        df.index = shiftdate(df.index, extra["tz"])
        del extra["tz"]

        # Fill extra columns with the same value everywhere
        # e.g., for coordinates when only specified in the header
        for e in extra:
            df[e] = extra[e]

        df["station"] = file_infos["stat"]
        df["network"] = file_infos["provider"]
        df["parameter"] = spec.lower()
        df["duration"] = default_duration
        df.rename(columns={spec.lower(): "obs"}, inplace=True)

        # If the parameter column is not a float, putting nans
        if df["obs"].dtype == "object":
            df["obs"] = np.nan

        # Removing rows with NaNs
        if kwargs.get("remove_nan", True):
            df = df[~np.isnan(df["obs"])]
            df = df[df["obs"] > na_values]
            df = df[df["obserror"] > err_na_values]

        # Rescales if needed
        if kwargs.get("rescale", False):
            coeffscale = rescale(obs_file, header)
            if np.isnan(coeffscale):
                info("Unknown scale, please check with provider")

            df["obs"] *= coeffscale
            df["obserror"] *= coeffscale

        # Converts unit
        df = convert_unit(df, ["obs", "obserror"], default_unit=default_unit)

        return df
