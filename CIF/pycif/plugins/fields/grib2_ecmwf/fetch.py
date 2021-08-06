import os

from pycif.utils import path
from .utils import find_valid_file


def fetch(ref_dir, ref_file, input_dates, target_dir, tracer=None, **kwargs):

    list_files = {}
    list_dates = {}
    for datei in input_dates:
        tmp_files = []
        tmp_dates = []
        for dd in input_dates[datei]:
            dir_dd = dd.strftime(ref_dir)
            files_3d, dates_3d = find_valid_file(dir_dd, tracer.file_3d, dd)
            tmp_files.extend(files_3d)
            tmp_dates.extend(dates_3d)

        # Fetching
        local_files = []
        for f in tmp_files:
            target_file = "{}/{}".format(target_dir, os.path.basename(f))
            path.link(f, target_file)
            local_files.append(target_file)

        list_files[datei] = list(set(local_files))
        list_dates[datei] = list(set(tmp_dates))

    return list_files, list_dates
