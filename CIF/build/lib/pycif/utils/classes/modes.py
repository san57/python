from types import MethodType

from pycif.utils.check.errclass import PluginError
from .baseclass import Plugin


class Mode(Plugin):
    def initiate_template(self):
        module = super(Mode, self).initiate(plg_type="mode")

        # Replacing auxiliary functions
        if hasattr(module, "execute"):
            self.execute = MethodType(module.execute, self)

    def execute(self, *args, **kwargs):
        """Default empty execute method"""

        raise PluginError("This is the default empty execute method")

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

        super(Mode, cls).register_plugin(
            name, version, module, plugin_type="mode"
        )
