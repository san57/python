from .fetch import fetch
from .read import read
from .write import write

requirements = {
    "domain": {"name": "CHIMERE", "version": "std", "empty": False}
}


def ini_data(plugin, **kwargs):

    # Default file names for CHIMERE: BOUN_CONCS
    if not hasattr(plugin, "file"):
        plugin.file = "BOUN_CONCS.%Y%m%d%H.{nho}.nc"
