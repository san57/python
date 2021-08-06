#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import absolute_import

from pycif.utils.classes.meteos import Meteo
from . import chimere_meteo
from . import dummy_csv
from . import lmdz_massflx

Meteo.register_plugin("LMDZ", "mass-fluxes", lmdz_massflx)
Meteo.register_plugin("dummy", "csv", dummy_csv)
Meteo.register_plugin("CHIMERE", "std", chimere_meteo)
