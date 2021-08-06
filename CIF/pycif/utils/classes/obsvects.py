from types import MethodType

import pandas as pd

from pycif.utils.check.errclass import PluginError
from .baseclass import Plugin


class ObsVect(Plugin):
    def __init__(self, **kwargs):
        """Create a Model Class"""

        super(ObsVect, self).__init__(**kwargs)

        # Attribute an empty datastore
        self.datastore = pd.DataFrame({})

    def initiate_template(self):
        module = super(ObsVect, self).initiate(plg_type="obsvect")

        # Replacing auxiliary functions:
        if hasattr(module, "init_y0"):
            self.init_y0 = MethodType(module.init_y0, self)

        if hasattr(module, "init_invprod"):
            self.init_invprod = MethodType(module.init_invprod, self)

        if hasattr(module, "rinvprod"):
            self.rinvprod = MethodType(module.rinvprod, self)

        if hasattr(module, "native2obsvect"):
            self.native2obsvect = MethodType(module.native2obsvect, self)

        if hasattr(module, "obsvect2native"):
            self.obsvect2native = MethodType(module.obsvect2native, self)

    def init_y0(self, *args, **kwargs):
        """Default empty init_y0 method"""

        raise PluginError("This is the default empty init_y0 method")

    def init_invprod(self, *args, **kwargs):
        """Default empty init_invprod method"""

        raise PluginError("This is the default empty init_invprod method")

    def rinvprod(self, *args, **kwargs):
        """Default empty rinvprod method"""

        raise PluginError("This is the default empty rinvprod method")

    def native2obsvect(self, *args, **kwargs):
        """Default empty native2obsvect method"""

        raise PluginError("This is the default empty native2obsvect method")

    def obsvect2native(self, *args, **kwargs):
        """Default empty obsvect2native method"""

        raise PluginError("This is the default empty obsvect2native method")

    @classmethod
    def register_plugin(cls, name, version, module, **kwargs):
        """Register a module for a plugin and version with possibly options

        Args:
            name (str):  name of the plugin
            version (str):  version of the plugin
            module (types.ModuleType): module defining the interface
                between pyCIF and the plugin
            plugin_type (str): type of plugin
            **kwargs (dictionary): default options for module

        """

        super(ObsVect, cls).register_plugin(
            name, version, module, plugin_type="obsvect"
        )
