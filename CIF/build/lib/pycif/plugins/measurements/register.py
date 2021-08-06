#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import absolute_import

from pycif.utils.classes.measurements import Measurement
from . import random
from . import standard

Measurement.register_plugin("standard", "std", standard)
Measurement.register_plugin("random", "std", random)
