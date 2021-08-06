from types import MethodType

from pycif.utils.check.errclass import PluginError
from .baseclass import Plugin


class Minimizer(Plugin):
    def initiate_template(self):
        module = super(Minimizer, self).initiate(plg_type="minimizer")

        # Replacing auxiliary functions
        if hasattr(module, "minimize"):
            self.minimize = MethodType(module.minimize, self)

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

        super(Minimizer, cls).register_plugin(
            name, version, module, plugin_type="minimizer"
        )

    def minimize(self, *args, **kwargs):
        """Default empty minimize function

        """
        raise PluginError("The function minimize was not defined")
