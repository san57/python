from .fetch import fetch
from .read import read
from .write import write

requirements = {
    "domain": {"name": "CHIMERE", "version": "std", "empty": False},
    "chemistry": {"name": "CHIMERE", "version": "gasJtab", "empty": False},
}


def ini_data(plugin, **kwargs):
    """Initializes the control vector from information in the Yaml file

    Args:
        plugin (pycif.classes.plugins): the plugin to initialize

    Return:
        - xb (explicitly and stored)
        - B (std and covariance definition), not stored
        - projectors and adjoints
        - product with B1/2


    """
    
    # Default type of emissions
    if not hasattr(plugin, "emis_type"):
        plugin.emis_type = "anthro"
    
    # Default file names for CHIMERE: AEMISSIONS
    if not hasattr(plugin, "file"):
        plugin.file = (
            "AEMISSIONS.%Y%m%d%H.{nho}.nc" if plugin.emis_type == "anthro"
            else "BEMISSIONS.%Y%m%d%H.{nho}.nc"
        )
