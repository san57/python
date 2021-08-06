################
Preparing inputs
################


:bash:`native2inputs`
---------------------

As mentioned earlier, we assume that all inputs and outputs are already computed here.
So, :bash:`native2inputs` only needs to copy or link the pre-computed file to the working directory.
:bash:`native2inputs` is based on the following structure:

.. code-block:: python

    def native2inputs(self, datastore, input_type, datei, datef, mode,
                      runsubdir, workdir):
        """Converts data at the model data resolution to model compatible input
        files.

        Args:
            self: the model Plugin
            input_type (str): type of input; should be consistent
                with self.required_inputs as specified in __init__.py
            datastore: data to convert
            datei, datef: date interval of the sub-simulation
            mode (str): running mode: one of 'fwd', 'adj' and 'tl'
            runsubdir (str): sub-directory for the current simulation
            workdir (str): the directory of the whole pyCIF simulation

        """

        if input_type == 'input_type_1':
            do_something()

        if input_type == 'input_type_2':
            do_something_else()

        ...

It is important to note that :bash:`native2inputs` needs that the :bash:`model`
has an attribute called :bash:`required_inputs`
defined in the :bash:`__init__.py` file:

.. code-block:: python

    def ini_data(plugin, **kwargs):
        """Initializes the model

        Args:
            plugin (dict): dictionary defining the plugin
            **kwargs (dictionary): possible extra parameters

        Returns:
            loaded plugin and directory with executable

        """

        info("Initializing the model")

        ...

        # Required inputs for running a simulation
        plugin.required_inputs = ['exe', 'nml', 'fluxes',
                                  'meteo', 'inicond', 'latcond', 'topcond',
                                  'any-input-type', 'your-output-for-the-tuto-purpose'...]

Once :bash:`native2inputs` is coded, it is necessary to test that it correctly sets up a working directory





