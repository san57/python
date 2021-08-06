##################################################
Pre-process observations to fit in model realm
##################################################

.. role:: bash(code)
   :language: bash


Observations cannot be compared to simulations directly as models simulate a discretized version of
the real world.
Prior to any simulation, pyCIF pre-computes which grid cells and time steps an observation data point should be compared to.

It consists in filling the following columns of the observation datastore:

    * i
    * j
    * level
    * tstep
    * tstep_glo
    * dtstep

Details on these columns are given :doc:`here<../../documentation/monitor>`.

This step is done automatically by the pycif class :bash:`obsvect`.
It is possible to do the processing using the Yaml configuration :doc:`here<Yaml_obsvect>`.
Try running :bash:`python -m pycif your_obsvect_yaml.yml`.

You should be returned the following error:

.. code-block:: python

    Exception: <pycif.utils.classes.obsvects.ObsVect object at 0x7f3e78d7b110>
    needs a Plugin 'model' to run properly
    there is none in its children nor at the level 0 of Yaml
    Please check your Yaml

This error is returned when a given :bash:`plugin` has :bash:`requirements` that are not fulfilled.
Further details :doc:`here<../../documentation/plugins/dependencies>`.

In the present case, the :bash:`obsvect` plugin requires a :bash:`model` (yours!) to work.
Therefore, you need to add a :bash:`model` paragraph to your Yaml configuration file:

.. code-block:: yaml

    # Transport model
    # Accepts any model registered in pycif.models
    model :
      plugin:
        name    : your-model
        version : your-model-version

      # Some other attributes

Examples of Yaml configuration attributes for :bash:`model` plugins are given
:doc:`here</configuration/models/index>`.


Temporal fit to the model time steps
------------------------------------

Some models can directly compare observations at any time with simulations.
However, most models compute simulations along time steps, and long simulations are often
split into a chain of sub-simulations.

The :bash:`obsvect` plugin automatically computes the time steps corresponding to a given observation
once the following attributes are associated to the :bash:`model` plugin:

- :bash:`subsimu_dates`:
    List of the starting dates of all the model sub-simulations chained to form the simulation window;
    it should include the end date of the chain simulation as well;
    for instance, LMDZ simulations are split into monthly simulations;
    :bash:`subsimu_dates` would then include the list of the first day of every months in the simulation window;
    if your model is not split into sub-simulations, :bash:`subsimu_dates` is a list with a single element,
    the starting date of the overall simulation window

- :bash:`tstep_all`:
    a list storing the starting dates of all the model time steps in the chain of simulation, including the end date of the chain simulation;
    for instance, LMDZ uses bi-hourly time steps; thus,  :bash:`tstep_all` will be a list with one value every 30 minutes
    in the simulation window

- :bash:`tstep_dates`:
    a dictionary with the elements of :bash:`subsimu_dates` as keys and the corresponding elements of :bash:`tstep_all`;
    for LMDZ, it will look like:

         .. code-block:: python

            tstep_dates = {datetime.datetime(2019, 1, 1):
                                [datetime.datetime(2019, 1, 1, 0, 0),
                                 datetime.datetime(2019, 1, 1, 0, 30),
                                 datetime.datetime(2019, 1, 1, 1, 0),
                                 ...
                                 datetime.datetime(2019, 1, 31, 23, 30),
                                 datetime.datetime(2019, 2, 1, 0, 0)],
                           datetime.datetime(2019, 2, 1):
                                [datetime.datetime(2019, 2, 1, 0, 0),
                                 datetime.datetime(2019, 2, 1, 0, 30),
                                 datetime.datetime(2019, 2, 1, 1, 0),
                                 ...
                                 datetime.datetime(2019, 2, 28, 23, 30),
                                 datetime.datetime(2019, 3, 1, 0, 0)],
                           ...}


To initialize these variables, it is necessary to attach a routine called :bash:`ini_periods.py` to your :bash:`model` plugin.
To do so, you need to code the routine in your :bash:`model` python module and in the corresponding :bash:`__init__.py` script,
add the line:

.. code-block:: python

    from ini_periods import ini_periods

With that being done, the :bash:`obsvect` plugin can use the above-mentioned variables as follows:

.. code-block:: python

    def tstep(self, *args, **kwargs):
        subsimu_dates = self.model.subsimu_dates

        ...

Horizontal fit to model domain
--------------------------------

pyCIF compares observations to simulations in corresponding grid cells.
To avoid recomputing over and over the corresponding grid cells, this information is calculated once for all during the
initialization of the :bash:`obsvect` plugin.
To allow the :bash:`obsvect` plugin to compute observation grid cells, the :bash:`model` plugin needs to be initialized
with a :bash:`domain` plugin attached to it.

To do this, please follow steps explained :doc:`here</devtutos/newplugin/newplugin>` and :doc:`here<../../documentation/plugins/dependencies>`.

.. hint::

    If this step is not applicable to your model (for instance, if the model compares simulation with the exact observation location),
    it is possible to use a virtual domain with only one grid cell:

    .. code-block:: yaml

        #####################################################################
        #####################################################################
        # Domain definition
        domain :
          plugin :
            name : dummy
            version : std
          xmin: -180
          xmax: 180
          nlon: 1
          ymin: -90
          ymax: 90
          nlat: 1
        #####################################################################
        #####################################################################



