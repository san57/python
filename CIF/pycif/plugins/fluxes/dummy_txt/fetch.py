import os

from pycif.utils import path
from .make import make


def fetch(ref_dir, ref_file, input_dates, target_dir, tracer=None, **kwargs):
    if ref_dir == "" and hasattr(tracer, "dir"):
        input_dir = tracer.dir
    else:
        input_dir = ref_dir

    if ref_file == "" and hasattr(tracer, "file"):
        input_file = tracer.file
    else:
        input_file = ref_file

    list_files = {
        datei: datei.strftime("{}/{}".format(input_dir, input_file))
        for datei in input_dates
    }
    list_dates = {datei: [] for datei in input_dates}

    # Fetching or creating if necessary
    for datei in input_dates:
        f = list_files[datei]
        target_file = "{}/{}".format(target_dir, os.path.basename(f))
        if os.path.isfile(f):
            path.link(f, target_file)

        else:
            flx = make(tracer, target_file, tracer.flx_text)
            tracer.write("", target_file, flx[0, 0])

        list_files[datei] = target_file

    return list_files, list_dates
