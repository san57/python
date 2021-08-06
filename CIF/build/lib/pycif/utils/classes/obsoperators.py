import datetime
from types import MethodType

from .baseclass import Plugin


class ObsOperator(Plugin):
    def obsoper(
        self,
        inputs,
        mode,
        run_id=0,
        datei=datetime.datetime(1979, 1, 1),
        datef=datetime.datetime(2100, 1, 1),
        workdir="./",
        **kwargs
    ):
        """The observation operator.
        This function maps information from the control space to the observation
        space and conversely depending on the running mode.

        Args:
            self (ObsOperator): the ObsOperator plugin
            inputs (Controlvect or Obsvect): the inputs of the fwd or adj mode
            mode (str): the running mode
            run_id (int): the ID of the current run (determines the
            sub-directory name
            datei (datetime.datetime): beginning of the simulation window
            datef (datetime.datetime): end of the simulation window
            workdir (str): path to the parent directory
            """
        return

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

        super(ObsOperator, cls).register_plugin(
            name, version, module, plugin_type="obsoperator"
        )

    def initiate_template(self):
        module = super(ObsOperator, self).initiate(plg_type="obsoperator")

        # Replacing auxiliary functions:
        if hasattr(module, "obsoper"):
            self.obsoper = MethodType(module.obsoper, self)
