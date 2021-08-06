import os

from pycif.utils import path
from pycif.utils.check import info


def default_fetch(
    ref_dir, ref_file, input_dates, target_dir, tracer=None, **kwargs
):

    # Picking info from tracer if not in inputs
    if ref_dir == "" and hasattr(tracer, "dir"):
        input_dir = tracer.dir
    else:
        input_dir = ref_dir

    if ref_file == "" and hasattr(tracer, "file"):
        input_file = tracer.file
    else:
        input_file = ref_file

    info("Fetching input files using directory and file format")
    info("{}/{}".format(input_dir, input_file))

    list_files = {}
    list_dates = {}
    for datei in input_dates:
        tmp_files = []
        tmp_dates = []
        for dd in input_dates[datei]:
            dir_dd = dd.strftime(input_dir)
            file_dd = dd.strftime(input_file)
            tmp_files.append("{}/{}".format(dir_dd, file_dd))

        # Fetching
        local_files = []
        for f in tmp_files:
            target_file = "{}/{}".format(target_dir, os.path.basename(f))
            path.link(f, target_file)
            local_files.append(target_file)

        list_files[datei] = list(set(local_files))
        list_dates[datei] = list(set(tmp_dates))

    return list_files, list_dates
