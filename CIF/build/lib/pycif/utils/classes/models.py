from types import MethodType

from pycif.utils.check.errclass import PluginError
from .baseclass import Plugin


class Model(Plugin):
    def __init__(self, **kwargs):
        """Create a Model Class"""

        super(Model, self).__init__(**kwargs)

        self.requirements.update({"domain": {"any": True, "empty": False}})

        # Default attributes
        self.required_inputs = []
        
        # Empty backup comps
        self.backup_comps = {}

    def initiate_template(self):
        module = super(Model, self).initiate(plg_type="model")

        # Replacing auxiliary functions
        if hasattr(module, "run"):
            self.run = MethodType(module.run, self)

        if hasattr(module, "native2inputs"):
            self.native2inputs = MethodType(module.native2inputs, self)

        if hasattr(module, "outputs2native"):
            self.outputs2native = MethodType(module.outputs2native, self)

        if hasattr(module, "make_input"):
            self.make_input = MethodType(module.make_input, self)

        if hasattr(module, "flushrun"):
            self.flushrun = MethodType(module.flushrun, self)

        if hasattr(module, "ini_periods"):
            self.ini_periods = MethodType(module.ini_periods, self)

        if hasattr(module, "compile"):
            self.compile = MethodType(module.compile, self)

        if hasattr(module, "ini_mapper"):
            self.ini_mapper = MethodType(module.ini_mapper, self)

        # Force initializing periods in ini_data
        if hasattr(self, "ini_data"):
            self.ini_data_orig = self.ini_data

            def ini_data(self, **kwargs):
                self.ini_data_orig(**kwargs)
                self.ini_periods(**kwargs)

            self.ini_data = MethodType(ini_data, self)

    def run(self, *args, **kwargs):
        """Default empty run method"""

        raise PluginError("This is the default empty run method")

    def native2inputs(self, *args, **kwargs):
        """Default empty native2inputs method"""

        raise PluginError("This is the default empty native2inputs method")

    def outputs2native(self, *args, **kwargs):
        """Default empty outputs2native method"""

        raise PluginError("This is the default empty outputs2native method")

    def make_input(self, *args, **kwargs):
        """Default empty make_input method"""

        raise PluginError("This is the default empty make_input method")

    def flushrun(self, *args, **kwargs):
        """Default empty flushrun method"""

        raise PluginError("This is the default empty flushrun method")

    def ini_periods(self, *args, **kwargs):
        """Default empty ini_periods method"""

        raise PluginError("This is the default empty ini_periods method")

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

        super(Model, cls).register_plugin(
            name, version, module, plugin_type="model"
        )

    def ini_mapper(self, *args, **kwargs):
        """Default empty ini_mapper method"""

        raise PluginError("This is the default empty ini_mapper method")
