#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import absolute_import

from pycif.utils.classes.fluxes import Fluxes
from . import chimere
from . import dummy_nc
from . import dummy_txt
from . import edgar_v5
from . import lmdz_bin
from . import lmdz_sflx
from . import flexpart

Fluxes.register_plugin("LMDZ", "sflx", lmdz_sflx)
Fluxes.register_plugin("LMDZ", "bin", lmdz_bin)
Fluxes.register_plugin("dummy", "nc", dummy_nc)
Fluxes.register_plugin("dummy", "txt", dummy_txt)
Fluxes.register_plugin("CHIMERE", "AEMISSIONS", chimere)
Fluxes.register_plugin("EDGAR", "v5", edgar_v5)
Fluxes.register_plugin("FLEXPART", "nc", flexpart)
