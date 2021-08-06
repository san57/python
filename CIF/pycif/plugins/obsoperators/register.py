#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pycif.utils.classes.obsoperators import ObsOperator
from . import standard

ObsOperator.register_plugin("standard", "std", standard)
