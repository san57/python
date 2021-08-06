import os
from pycif.utils import path


def fetch(ref_dir, ref_file, input_dates, target_dir, tracer=None, **kwargs):
    if ref_dir == '' and hasattr(tracer, 'dir'):
        input_dir = tracer.dir
    else:
        input_dir = ref_dir

    if ref_file == '' and hasattr(tracer, 'file'):
        input_file = tracer.file
    else:
        input_file = ref_file

    list_files = {datei: datei.strftime('{}/{}'.format(input_dir, input_file))
                  for datei in input_dates}
    
    if hasattr(tracer, "file_glob"):
        list_files_glo = \
            {datei: datei.strftime('{}/{}'.format(input_dir, tracer.file_glob))
             for datei in input_dates}

    list_dates = {datei: input_dates[datei] for datei in input_dates}
    
    # Fetching
    for datei in input_dates:
        f = list_files[datei]
        target_file = '{}/{}'.format(target_dir, os.path.basename(f))
        path.link(f, target_file)
        list_files[datei] = target_file
        
        # Initialize global file for later
        if hasattr(tracer, "file_glob"):
            f = list_files_glo[datei]
            target_file = '{}/{}'.format(target_dir, os.path.basename(f))
            path.link(f, target_file)
    
    return list_files, list_dates
