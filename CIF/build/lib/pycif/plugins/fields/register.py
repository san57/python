#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import absolute_import

from pycif.utils.classes.fields import Fields
from . import chimere_icbc
from . import grib2_ecmwf
from . import lmdz_ic
from . import lmdz_prescrconcs

Fields.register_plugin("CHIMERE", "icbc", chimere_icbc)
Fields.register_plugin("LMDZ", "ic", lmdz_ic)
Fields.register_plugin("LMDZ", "prescrconcs", lmdz_prescrconcs)
Fields.register_plugin("ECMWF", "grib2", grib2_ecmwf)
