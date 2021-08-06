#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os

from pycif.utils import path


def make_meteo(self, data, ddi, ddf, runsubdir, mode):
    datei = min(ddi, ddf)

    dir_meteo = data.datastore[("meteo", "")]["dirorig"]
    # Links meteorological mass fluxes
    for ftype in ["defstoke", "fluxstoke", "fluxstokev", "phystoke"]:
        meteo_file = "{}.an{}.m{:02d}.nc".format(
            ftype, datei.year, datei.month
        )
        source = "{}/{}".format(dir_meteo, meteo_file)
        target = "{}/{}.nc".format(runsubdir, ftype)

        if ftype == "defstoke" and not os.path.isfile(source):
            meteo_file = "{}.nc".format(ftype)
            source = "{}/{}".format(dir_meteo, meteo_file)

        path.link(source, target)

    # Links vertical coordinates from previous simulations
    source = self.file_vcoord
    target = "{}/vcoord.nc".format(runsubdir)
    path.link(source, target)
