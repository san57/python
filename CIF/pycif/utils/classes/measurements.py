from types import MethodType

from pycif.utils.check.errclass import PluginError
from pycif.utils.datastores.empty import init_empty
from .baseclass import Plugin


class Measurement(Plugin):
    def initiate_template(self):
        module = super(Measurement, self).initiate(plg_type="measurements")

        # Replacing auxiliary functions
        if hasattr(module, "parse_tracers"):
            self.parse_tracers = MethodType(module.parse_tracers, self)

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

        super(Measurement, cls).register_plugin(
            name, version, module, plugin_type="measurements"
        )

    def ini_data(self, **kwargs):
        """Initializes the measurement plugin. This is pased on iterating over
        multiple species and providers

        Args:
            plugin (MeasurementPlugin): Measurement definition

        Returns:
            Updates on the fly the measurements
        """

        datei = self.datei
        datef = self.datef

        # If the measurement definition is empty in the Yaml,
        # return an empty datastore
        if not hasattr(self, "species"):
            self.datastore = init_empty()

        else:
            self.datastore = self.parse_tracers(datei, datef, **kwargs)

    def parse_tracers(self, *args, **kwargs):
        """Default empty minimize function

        """
        raise PluginError("The function parse_tracers was not defined")
