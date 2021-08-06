import datetime
import glob
import os

import numpy as np

from pycif.utils import path
from pycif.utils.check import info


def fetch(ref_dir, ref_file, input_dates, target_dir, tracer=None, **kwargs):

    list_files = {}
    list_dates = {}
    for datei in input_dates:
        tmp_files = []
        tmp_dates = []
        for dd in input_dates[datei]:
            file_flx = dd.strftime(ref_file)
            dir_flx = dd.strftime(ref_dir)
            date_flx = dd

            if not os.path.isfile(
                "{}/{}".format(dir_flx, file_flx)
            ) and getattr(tracer, "closest_year", False):
                info(
                    "Warning: could not find correct year for EDGAR; "
                    "using closest available one"
                )
                list_dates_avail = [
                    datetime.datetime.strptime(os.path.basename(f), ref_file)
                    for f in glob.glob("{}/v50_*nc".format(dir_flx))
                ]
                delta_dates = np.abs(dd - np.array(list_dates_avail))
                date_flx = list_dates_avail[np.argmin(delta_dates)]
                file_flx = date_flx.strftime(ref_file)

            tmp_files.append("{}/{}".format(dir_flx, file_flx))
            tmp_dates.append(date_flx)

        # Fetching
        local_files = []
        for f in tmp_files:
            target_file = "{}/{}".format(target_dir, os.path.basename(f))
            path.link(f, target_file)
            local_files.append(target_file)

        list_files[datei] = list(set(local_files))
        list_dates[datei] = list(set(tmp_dates))

    return list_files, list_dates
