#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import absolute_import

from pycif.utils.classes.obsparsers import ObsParser
from . import wdcgg

ObsParser.register_parser("WDCGG", "std", wdcgg)
