from types import MethodType

from pycif.utils.check.errclass import PluginError
from .baseclass import Plugin


class Transform(Plugin):
    def initiate_template(self):
        module = super(Transform, self).initiate(plg_type="transform")

        # Replacing auxiliary functions:
        if hasattr(module, "ini_mapper"):
            self.ini_mapper = MethodType(module.ini_mapper, self)
        #
        # if hasattr(module, 'native2state'):
        #     self.native2state = MethodType(module.native2state, self)
        #
        # if hasattr(module, 'obsvect2native'):
        #     self.obsvect2native = MethodType(module.obsvect2native, self)
        #
        # if hasattr(module, 'native2obsvect'):
        #     self.native2obsvect = MethodType(module.native2obsvect, self)

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

        super(Transform, cls).register_plugin(
            name, version, module, plugin_type="transform"
        )

    @classmethod
    def get_transform(cls, plg):
        """Get the correct Parser for a provider and file_format_id

        Args:
            provider (str):  provider of the input file
            file_format_id (str): name of the type of file with a given format

        Returns:
            Parser: Parser for provider and file_format_id
        """

        return cls.load_registered(
            plg.orig_name, "std", "transform", plg_orig=plg
        )

    def ini_mapper(self, *args, **kwargs):
        """Default empty ini_mapper method"""

        raise PluginError("This is the default empty ini_mapper method")
