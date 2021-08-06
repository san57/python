#################################
Comparing outputs to observations
#################################

.. role:: bash(code)
   :language: bash

Once observations are implemented in the model, the next step is to compare them with model outputs.
The objective of the present step is to "run" pyCIF in :doc:`forward mode</configuration/modes/forward>`
and produce a file comparing prescribed observations with simulations.
It is assumed that all inputs and outputs from a pre-existing simulations are available for your model.

A :bash:`forward` mode requires the following :bash:`model` methods to be provided,
as they are called by the :bash:`obsoper` plugin:

1. :bash:`native2inputs`:
    Prepare all inputs for the sub-simulation of interest

2. :bash:`outputs2native`:
    Read model outputs into pyCIF variables

3. :bash:`run`:
    Run the model

---------------------
:bash:`native2inputs`
---------------------

As mentioned earlier, we assume that all inputs and outputs are already computed.
So, :bash:`native2inputs` does not need to do anything.
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

        # do nothing
        return


---------------------
:bash:`run`
---------------------

We do not try to really run the model at this step.
Therefore, :bash:`run` will only copy pre-computed outputs to the current working directory.

.. code-block:: python

    def run(self, runsubdir, mode, workdir, do_simu=True, **kwargs):
        """Run the model

        Args:
            self: the model Plugin
            runsubdir (str): working directory for the current run
            mode (str): forward or backward
            workdir (str): pyCIF working directory
            do_simu (bool): re-run or not existing simulation

        """

        # Copy or link your outputs to runsubdir


----------------------
:bash:`outputs2native`
----------------------

At this step, raw model outputs are processed and put into pyCIF variables.
For a forward run, the objective is to recompose observation equivalents.
pyCIF expects simulations to be included into a pandas DataFrame.

.. code-block:: python


    def outputs2native(self, runsubdir, mode='fwd', dump=True):
        """Reads outputs to pyCIF objects.

        If the mode is 'fwd' or 'tl', only onservation-like outputs are extracted.
        For the 'adj' mode, all outputs relative to model sensitivity are extracted.

        Dumps to a NetCDF file with output concentrations if needed"""

        if not hasattr(self, 'dataobs'):
            self.dataobs = init_empty()

        if not hasattr(self, 'datasensit'):
            self.datasensit = {}

        if mode in ['fwd', 'tl']:
            # If no simulated concentration is available just pass
            file = "{}/mod.txt".format(runsubdir)
            if os.stat(file).st_size == 0:
                info("CHIMERE ran without any observation to be compared with for sub-simu "
                "only CHIMERE's outputs are available")
                self.dataobs.loc[:, 'sim'] = np.nan
                return

            # Read simulated concentrations
            sim = pd.read_csv(file, sep='\s+', header=None)[6]

            # Loop over observations in active species
            mask = self.dataobs['parameter'].str.upper() \
                .isin(list(self.chemistry.species['name']))

            # Putting values to the local data store
            # Assumes arithmetic averages upto now
            inds = [0] + list(np.cumsum(self.dataobs.loc[mask, 'dtstep'][:-1]))

            column = 'sim' if mode == 'fwd' else 'sim_tl'
            self.dataobs.loc[mask, column] = \
                [sim.iloc[k:k + dt].sum()
                 for k, dt in zip(inds, self.dataobs.loc[mask, 'dtstep'])]

            return self.dataobs















