# -*- coding: utf-8 -*-

from __future__ import absolute_import

import os

import numpy as np
from dateutil.tz import tzoffset


def remap_head(s):
    """Adapts names to extract values in WDCGG header

    """

    if s.lower() == "lat":
        return "latitude:"

    elif s.lower() == "lon":
        return "longitude:"

    elif s.lower() == "alt":
        return "altitude:"

    elif s.lower() == "unit":
        return "measurement unit"

    elif s.lower() == "tz":
        return "time zone"

    else:
        return s.lower()


def find_header(id_extract, header):
    """Finds the value of a constant parameter (e.g., latitude, altitude, etc.)
    in the header of a file

    """

    for ln in header:
        if remap_head(id_extract) in ln.lower():
            value = ln.lower().split(":")[1].strip()
            try:
                return float(value)
            except ValueError:
                return value

    raise ValueError("Couldn't extract {}".format(id_extract))


def rescale(obs_file, header):
    """Finds out on what scale the measurement was reported and returns the
    corresponding scaling factor.

    Notes:
        If the scale is not in the list of recognized scales, then returns a
        NaN to kill the data

    """

    try:
        scale = find_header("scale", header)

    except BaseException:
        scale = parse_file(obs_file)["provider"]

    if "04" in scale or "wmo" in scale.lower():
        return 1.0

    elif scale == "CSIRO94":
        return 1.01219

    elif "NIST" in scale or "USA" in scale:
        return 0.998

    elif "tohoku" in scale.lower():
        return 1.0003

    elif "aircore" in scale.lower():
        return 1.0124

    elif scale == "Manufacture's":
        return 0.997

    elif "NIES" in scale:
        return 0.997

    # If the scale is not known, then returns NaN
    return np.nan


def parse_file_name(obs_file, **kwargs):
    """Parses WDCGG file name and extract corresponding information.

    This is based on WDCGG standard naming format as detailed in:
    http://ds.data.jma.go.jp/gmd/wdcgg/pub/data/WDCGG_filename_format.pdf

    """
    filesplit = os.path.basename(obs_file).split(".")

    infos = {}
    infos["stat"] = filesplit[0][:3]
    infos["provider"] = filesplit[1].replace("_", "-")
    infos["site category"] = filesplit[2]
    infos["obs type"] = filesplit[-5]
    infos["parameter"] = filesplit[-4]
    infos["freq"] = filesplit[-3]

    return infos


def convert_unit(df, params, unit="ppm", default_unit="ppm"):
    """Converts between ppb, ppm, ppt. Default is conversion to ppm

    """

    if "unit" in df.columns:
        for p in params:
            # Change missing unit to default unit
            df.loc[df["unit"] == "", "unit"] = default_unit

            # First conversion to ppm as a common reference unit
            df.loc[df["unit"] == "ppt", p] /= 1e6
            df.loc[df["unit"] == "ppb", p] /= 1e3
            df.loc[df["unit"] == "ppbv", p] /= 1e3
            df.loc[df["unit"] == "nmol.mol-1", p] /= 1e3
            df.loc[df["unit"] == "nmol.mol-¹", p] /= 1e3

            # Then conversion to target unit if needed
            if unit in ["ppb", "ppbv", "nmol.mol-1", "nmol.mol-¹"]:
                df[p] *= 1e3
            elif unit == "ppt":
                df[p] *= 1e6
            elif unit == "ppm":
                pass
            else:
                raise ValueError(unit + " is not a valid unit for conversion")

    df["unit"] = unit

    return df


def shiftdate(dates, tz):
    """Shits dates according to a time zone as define in WDCCGG files

    """

    if tz in ["utc", "utc+0", "utca+0", "utc +0", "", None]:
        return dates

    utc_code = [w for w in tz.split() if "utc" in w][0]

    shift = int(utc_code[4:]) * (-1 + 2 * (utc_code[3] == "+"))

    tzlocal = tzoffset("local", 60 * 60 * shift)
    dates = dates.tz_localize(tzlocal)
    dates = dates.tz_convert("UTC")
    dates = dates.tz_localize(None)

    return dates


def remap_extract(s):
    """Adapts names to extract columns in the context of WDCGG

    """

    if s.lower() == "obserror":
        return "sd"

    elif s.lower() == "flag":
        return "f"

    elif s.lower() == "mcf":
        return "ch3ccl3"

    else:
        return s.lower()
