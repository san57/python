#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pycif.utils.classes.domains import Domain
from . import chimere
from . import dummy
from . import lmdz
from . import flexpart

Domain.register_plugin("LMDZ", "std", lmdz)
Domain.register_plugin("CHIMERE", "std", chimere)
Domain.register_plugin("dummy", "std", dummy)
Domain.register_plugin("FLEXPART", "std", flexpart)

del lmdz
del chimere
