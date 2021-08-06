from .baseclass import Plugin


class Platform(Plugin):
    def initiate_template(self):
        module = super(Platform, self).initiate(plg_type="platform")

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

        super(Platform, cls).register_plugin(
            name, version, module, plugin_type="platform"
        )
