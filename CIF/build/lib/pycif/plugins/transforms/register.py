#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pycif.utils.classes.transforms import Transform
from . import (
    families,
    isotopes,
    satellites,
    timeavg,
    regrid,
    unit_conversion,
    time_interpolation,
    vertical_interpolation,
    run_model,
    init_parameter,
    dump_parameter,
)

Transform.register_plugin("families", "std", families)
Transform.register_plugin("isotopes", "std", isotopes)
Transform.register_plugin("satellites", "std", satellites)
Transform.register_plugin("timeavg", "std", timeavg)
Transform.register_plugin("regrid", "std", regrid)
Transform.register_plugin("unit_conversion", "std", unit_conversion)
Transform.register_plugin("time_interpolation", "std", time_interpolation)
Transform.register_plugin(
    "vertical_interpolation", "std", vertical_interpolation
)
Transform.register_plugin("run_model", "std", run_model)
Transform.register_plugin("init_parameter", "std", init_parameter)
Transform.register_plugin("dump_parameter", "std", dump_parameter)
