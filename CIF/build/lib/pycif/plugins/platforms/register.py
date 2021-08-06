#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import absolute_import

from pycif.utils.classes.platforms import Platform
from . import lsce_obelix
from . import tgcc_ccrt

Platform.register_plugin("LSCE", "obelix", lsce_obelix)
Platform.register_plugin("TGCC-CCRT", "std", tgcc_ccrt)
