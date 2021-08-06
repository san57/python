from types import MethodType

from pycif.utils.check.errclass import PluginError
from .baseclass import Plugin


class Chemistry(Plugin):
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

        super(Chemistry, cls).register_plugin(
            name, version, module, plugin_type="chemistry"
        )

    def read_chemicalscheme(self, *args, **kwargs):
        """Read a chemical scheme from an existing file

        Args:
            self (Chemistry): plugin defining the chemistry. Should include
            dirscheme grid to be able to read the chemistry from files

        Return:
            Characteristics of the chemical scheme
        """
        raise PluginError("The function read_scheme was not defined")

    def create_chemicalscheme(self, *args, **kwargs):
        """Creates a chemical scheme if needed

        Args:
            chemistry (dictionary): dictionary defining the chelical scheme
        """
        raise PluginError("The function create_chemicalscheme was not defined")

    def initiate_template(self, plg_type=None):
        """Initializes a Plugin given the template

        Args:
            plg_type (str): the plugin type

        """

        # Initializes the Plugin from the parent method
        module = super(Chemistry, self).initiate(plg_type="chemistry")

        # Replacing auxiliary functions
        if hasattr(module, "read_chemicalscheme"):
            self.read_chemicalscheme = MethodType(
                module.read_chemicalscheme, self
            )

        if hasattr(module, "create_chemicalscheme"):
            self.create_chemicalscheme = MethodType(
                module.create_chemicalscheme, self
            )
