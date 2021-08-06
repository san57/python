from types import MethodType

from .baseclass import Plugin


class Meteo(Plugin):
    def initiate_template(self):
        module = super(Meteo, self).initiate(plg_type="meteo")

        # Replacing auxiliary functions:
        if hasattr(module, "read"):
            self.read = MethodType(module.read, self)
        if hasattr(module, "write"):
            self.write = MethodType(module.write, self)
        if hasattr(module, "fetch") and not hasattr(self, "fetch"):
            self.fetch = MethodType(module.fetch, self)

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

        super(Meteo, cls).register_plugin(
            name, version, module, plugin_type="meteo"
        )
