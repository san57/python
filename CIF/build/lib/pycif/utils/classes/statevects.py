from types import MethodType

from pycif.utils.check.errclass import PluginError
from .baseclass import Plugin


class StateVect(Plugin):
    def initiate_template(self):
        module = super(StateVect, self).initiate(plg_type="statevect")

        # Replacing auxiliary functions:
        # Need some function to convert information from the state space
        # to the space of model inputs
        if hasattr(module, "state2inputs"):
            self.state2inputs = MethodType(module.state2inputs, self)

        if hasattr(module, "outputs2state"):
            self.outputs2state = MethodType(module.outputs2state, self)

        if hasattr(module, "control2native"):
            self.control2native = MethodType(module.control2native, self)

        if hasattr(module, "native2control"):
            self.native2control = MethodType(module.native2control, self)

        # Define how to multiple by sqrt-B, i.e. translating information
        # from the minimization space to the physical space
        if hasattr(module, "sqrtbprod"):
            self.sqrtbprod = MethodType(module.sqrtbprod, self)

        if hasattr(module, "sqrtbprod_ad"):
            self.sqrtbprod_ad = MethodType(module.sqrtbprod_ad, self)

        # Load and dump control vectors
        if hasattr(module, "dump"):
            self.dump = MethodType(module.dump, self)

        if hasattr(module, "load"):
            self.load = MethodType(module.load, self)

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

        super(StateVect, cls).register_plugin(
            name, version, module, plugin_type="statevect"
        )

    def state2inputs(self, *args, **kwargs):
        """Default empty state2inputs method"""

        raise PluginError("This is the default empty state2inputs method")

    def outputs2state(self, *args, **kwargs):
        """Default empty outputs2state method"""

        raise PluginError("This is the default empty outputs2state method")

    def control2native(self, *args, **kwargs):
        """Default empty control2native method"""

        raise PluginError("This is the default empty control2native method")

    def native2control(self, *args, **kwargs):
        """Default empty native2control method"""

        raise PluginError("This is the default empty native2control method")

    def sqrtbprod(self, *args, **kwargs):
        """Default empty sqrtbprod method"""

        raise PluginError("This is the default empty sqrtbprod method")

    def sqrtbprod_ad(self, *args, **kwargs):
        """Default empty sqrtbprod_ad method"""

        raise PluginError("This is the default empty sqrtbprod_ad method")

    def dump(self, *args, **kwargs):
        """Default empty dump method"""

        raise PluginError("This is the default empty dump method")

    def load(self, *args, **kwargs):
        """Default empty load method"""

        raise PluginError("This is the default empty load method")
