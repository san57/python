from __future__ import absolute_import

import glob
from os.path import dirname, basename, isdir

from .plugins import chemistries, domains, fields, fluxes, measurements, \
    meteos, minimizers, models, modes, obsoperators, obsparsers, obsvects, \
    platforms, simulators, statevect, transforms

modules = glob.glob(dirname(__file__) + "/*")

__all__ = [basename(f) for f in modules if isdir(f) and not f.endswith(".py")]
