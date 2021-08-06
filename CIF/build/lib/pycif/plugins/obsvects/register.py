#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pycif.utils.classes.obsvects import ObsVect
from . import standard

ObsVect.register_plugin("standard", "std", standard)
