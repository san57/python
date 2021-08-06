#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import absolute_import

from pycif.utils.classes.models import Model
from . import chimere
from . import dummy
from . import lmdz
from . import flexpart

Model.register_plugin("LMDZ", "std", lmdz)
Model.register_plugin("CHIMERE", "std", chimere)
Model.register_plugin("dummy", "std", dummy)
Model.register_plugin("FLEXPART", "std", flexpart)
