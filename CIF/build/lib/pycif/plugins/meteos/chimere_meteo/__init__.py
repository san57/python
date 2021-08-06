from __future__ import absolute_import

from .fetch import fetch
from .read import read

requirements = {
    "domain": {"name": "CHIMERE", "version": "std", "empty": False}
}


def ini_data(plugin, **kwargs):

    # Default file names for CHIMERE: METEO
    if not hasattr(plugin, "file"):
        plugin.file = "METEO.%Y%m%d%H.{nho}.nc"
