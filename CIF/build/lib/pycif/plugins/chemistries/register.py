#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pycif.utils.classes.chemistries import Chemistry
from . import chimere

Chemistry.register_plugin("CHIMERE", "gasJtab", chimere)

del chimere
