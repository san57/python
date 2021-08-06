############
Observations
############

.. role:: bash(code)
   :language: bash


pyCIF uses observations to optimize fluxes or other control variables,
as well as to simply compare to forward simulations.
Observations are stored in a Pandas DataFrame when loaded in pyCIF;
pyCIF can read or dump observations from/to compatible NetCDF files as described
below.

Observation files
^^^^^^^^^^^^^^^^^

pyCIF uses storage NetCDF :bash:`monitor.nc` files for observations inspired from the
`ObsPACK <https://www.esrl.noaa.gov/gmd/ccgg/obspack/>`__ standard format.

Below is the list of variables required by pyCIF:

:index:
    ID number of each observation point
:date:
    starting date
:duration:
    duration in hours
:station:
    station ID (max 3 characters)
:network:
    network name (max 7 characters)
:parameter:
    name of the observed parameter or species
:lon:
    longitude of the measurement
:lat:
    latitude of the measurement
:alt:
    altitude (in m a.s.l)
:i:
    latitudinal index of the grid point corresponding to lon/lat
:j:
    longitudinal index of the grid point corresponding to lon/lat
:level:
    level number in the model
:tstep:
    time-step number in the sub-simulation of the measurement
:tstep_glo:
    time-step number in the complete chain of model sub-simulations
:dtstep:
    number of time steps in the model over which the measurement spans
:obs:
    observed value
:obserror:
    error on the observation (used to defined the matrix :math:`\mathbf{R}`)
:sim:
    simulated value
:sim_tl:
    simulated increments in the tangent-linear model

Utility functions are available in pyCIF to read and dump observation files.


Observation  data structure in pyCIF
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

pyCIF handles observations as `pandas.DataFrames <https://pandas.pydata.org/pandas-docs/stable/reference/api/pandas.DataFrame.html>`__.
The python module pandas offers very powerful data manipulations.
One can find some details on pandas as used in the CIF :doc:`here<pandas>`

It is possible to import/export observations from/to a :bash:`monitor.nc` file with the following commands:

.. code:: python

        from pycif.utils.datastores import dump

        # Reading the monitor file into a Pandas dataframe
        datastore = dump.read_datastore(your_monitor_file)

        print datastore

        # Exporting to a NetCDF
        dump.dump_datastore(datastore, file_monit=your_monitor_file)

If one needs to generate manually a pyCIF datastore, it is possible to fill it manually:

.. code:: python

        from pycif.utils.datastores import empty
        from pycif.utils.datastores.dump import dump_datastore

        # Initializes an empty datastore
        datastore = empty.init_empty()

        # Fill the datastore manually

        # Save the datastore to some netCDF
        dump_datastore(datastore, file_monitor=/some/path/where/to/save)
