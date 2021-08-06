#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

from pycif.utils import path
from pycif.utils.check import info


def fetch(
    ref_dir,
    ref_file,
    input_dates,
    target_dir,
    tracer=None,
    filetypes=["defstoke", "fluxstoke", "fluxstokev", "phystoke"],
    **kwargs
):
    """Reads meteorology and links to the working directory

    Args:
        meteo (dictionary): dictionary defining the domain. Should include
        dirmeteo to be able to read the meteorology
        datei (datetime.datetime): initial date for the inversion window
        datef (datetime.datetime): end date for the inversion window
        workdir (str): path to the working directory where meteo files
                       should be copied
        logfile (str): path to the log file
        filetypes ([str]): list of file radicals to copy in the working
                           directory
        **kwargs (dictionary): extra arguments

    Return:
        ????????

    Notes: At some point, include option to compute mass fluxes for LMDz,
    with different physics What is needed to do that? Possible only on CCRT?
    Flexibility to define new domains Can be very heavy and not necessarily
    relevant

    """

    info("Copying meteo files from {} to {}".format(ref_dir, target_dir))

    # Create the sub-directory to store meteo files
    path.init_dir(target_dir)

    # Loop over dates and file types
    for date in input_dates:
        for filetype in filetypes:
            meteo_file = "{}.an{}.m{:02d}.nc".format(
                filetype, date.year, date.month
            )

            if filetype == "defstoke" and not os.path.isfile(
                ref_dir + meteo_file
            ):
                meteo_file = filetype + ".nc"

            target = "{}/{}".format(target_dir, meteo_file)
            source = "{}/{}".format(ref_dir, meteo_file)
            path.link(source, target)

    list_files = {datei: [] for datei in input_dates}
    list_dates = {datei: [] for datei in input_dates}

    return list_files, list_dates
