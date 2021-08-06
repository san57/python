#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os

from pycif.utils.netcdf import readnc


def read(
    self,
    name,
    tracdir,
    tracfile,
    varnames,
    dates,
    interpol_flx=False,
    tracer=None,
    model=None,
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

    for date in dates:
        for filetype in filetypes:
            meteo_file = "{}.an{}.m{:02d}.nc".format(
                filetype, date.year, date.month
            )

            if filetype == "defstoke" and not os.path.isfile(
                "{}/{}".format(tracdir, meteo_file)
            ):
                meteo_file = filetype + ".nc"

            target = "{}/{}".format(tracdir, meteo_file)

            # Loading information on time steps
            if filetype == "defstoke" and not hasattr(self, "offtstep"):
                vars = readnc(target, ["dtvr", "istdyn"])
                offtstep = vars[0][0, 0] * vars[1][0, 0]
                self.offtstep = offtstep
