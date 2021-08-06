import datetime
import glob
import os

#import cfgrib
import numpy as np

from pycif.utils.check import info


def grib_file_reader(filepath, varname, attribute=None):
    """
        Filepath is the absolute file path
        Dimension_i is the name of dimension i, e.i. given between cotes ''
        if there is not third dimension, Dimension_3's value = None
        Variable is the variable's name, e.i. given between cotes ''
    """

    # Forcing import of attributes if not in attribute key list
    if (
        attribute is not None
        and attribute not in cfgrib.dataset.DATA_ATTRIBUTES_KEYS
    ):
        cfgrib.dataset.DATA_ATTRIBUTES_KEYS.append(attribute)
        cfgrib.dataset.ALL_KEYS = sorted(cfgrib.dataset.ALL_KEYS + [attribute])

    if not os.path.exists(filepath):
        info("{} was not found".format(filepath))
        raise IOError

    info("Reading {}".format(filepath))

    df = cfgrib.open_file(filepath)

    if len(np.shape(varname)) == 0:
        varnames = [varname]
    else:
        varnames = varname[:]

    varout = []
    for name in varnames:
        var = df.variables[name].data

        if hasattr(var, "build_array"):
            var = var.build_array()
        varout.append(var)

    # Fetching attributes if needed
    if attribute is not None:
        for name in df.variables:
            attr = df.variables[name].attributes.get(
                "GRIB_{}".format(attribute), None
            )
            if attr is not None:
                return attr
        raise Exception(
            "Could not find attribute {} in {}".format(attribute, filepath)
        )

    if len(np.shape(varname)) == 0:
        return varout[0]
    else:
        return varout


def find_valid_file(ref_dir, file_format, dd):
    ref_pref = file_format.split(".")[0]
    list_files = np.array(glob.glob("{}/{}.*grb2".format(ref_dir, ref_pref)))

    ref_dates = np.array(
        [
            datetime.datetime.strptime(
                os.path.basename(f).split(".")[-2][:8], "%Y%m%d"
            )
            for f in list_files
        ]
    )
    ref_deltas = np.array(
        [int(os.path.basename(f).split(".")[-2][12:]) for f in list_files]
    )

    file_dates = np.array(
        [
            d + datetime.timedelta(hours=dlt)
            for d, dlt in zip(ref_dates, ref_deltas)
        ]
    )

    mask = (file_dates - dd) <= datetime.timedelta(0)
    file_ref1 = list_files[mask][np.argmax(file_dates[mask])]
    date_ref1 = file_dates[mask].max()

    mask = (file_dates - dd) >= datetime.timedelta(0)
    file_ref2 = list_files[mask][np.argmin(file_dates[mask])]
    date_ref2 = file_dates[mask].min()

    return ([file_ref1, file_ref2], [date_ref1, date_ref2])
