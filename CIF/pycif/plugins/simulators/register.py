#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pycif.utils.classes.simulators import Simulator
from . import dummy
from . import gausscost

Simulator.register_plugin("gausscost", "std", gausscost)
Simulator.register_plugin("dummy_txt", "std", dummy)
