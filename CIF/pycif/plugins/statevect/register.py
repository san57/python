#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pycif.utils.classes.statevects import StateVect
from . import standard

StateVect.register_plugin("standard", "std", standard)
