# -*- coding: utf-8 -*-
"""Contains all recognized state vectors.
Automatically loads all pre-defined models as sub-modules of pycif.models
"""
from __future__ import absolute_import

import glob
from os.path import dirname, basename, isdir

from . import register

modules = glob.glob(dirname(__file__) + "/*")

__all__ = [basename(f) for f in modules if isdir(f) and not f.endswith(".py")]
