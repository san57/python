from __future__ import absolute_import
import os
from types import MethodType

from pycif.utils import path
from pycif.utils.check import verbose

from .read import read
from .read_glob import read_glob
from .write import write
from .fetch import fetch

requirements = {'domain': {'name': 'FLEXPART', 'version': 'std', 'empty': False}}

