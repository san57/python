#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pycif.utils.classes.modes import Mode
from . import variational, analytic, forward, footprint, adjtl_test, postproc

Mode.register_plugin("4dvar", "std", variational)
Mode.register_plugin("analytic", "std", analytic)
Mode.register_plugin("forward", "std", forward)
Mode.register_plugin("footprint", "std", footprint)
Mode.register_plugin("adj-tl_test", "std", adjtl_test)
Mode.register_plugin("post-proc", "std", postproc)
