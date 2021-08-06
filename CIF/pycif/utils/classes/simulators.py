from types import MethodType

from pycif.utils.check.errclass import PluginError
from .baseclass import Plugin


class Simulator(Plugin):
    def initiate_template(self):
        module = super(Simulator, self).initiate(plg_type="simulator")

        # Replacing auxiliary functions:
        if hasattr(module, "simul"):
            self.simul = MethodType(module.simul, self)

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

        super(Simulator, cls).register_plugin(
            name, version, module, plugin_type="simulator"
        )

    def simul(self, *args, **kwargs):
        """Default empty simul method"""

        raise PluginError("This is the default empty simul method")
