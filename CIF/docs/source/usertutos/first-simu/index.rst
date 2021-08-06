#####################################
How to run a first forward simulation
#####################################

.. role:: bash(code)
   :language: bash


Here we generate pyCIF outputs for a forward run using
:doc:`the Toy Gaussian Model</documentation/plugins/models/toy-gaussian>`.
To do so, follow the steps.
Please note :doc:`here</documentation/paths>` some instructions about how pyCIF deals with paths.

1. Set up a Yaml configuration:

    - copy the content of the reference Yaml file available :doc:`here<yaml_fwd>` to a directory of your choice
    - modify paths to fit your own configuration:
        - l. 31: :bash:`workdir`: where outputs will be generated
        - l. 76: :bash:`file_pg`: Pasquil-Gifford data, provided with CIF codes
        - l. 98: :bash:`file_monitor`: randomly generated observations
        - l. 135: :bash:`file_obsvect`: observations fitted to the model grid; does not need to be the same as :bash:`file_monitor`
        - l. 184-189: :bash:`dir`: directory where fluxes and correlations will be generated
        - l. 199: :bash:`dir`: directory to example meteo files, provided with CIF codes; if does not exist, will generate random meteo parameters

2. Run pyCIF based on your Yaml file:

    .. code-block:: bash

        python -m pycif path-to-your-yaml

3. Explore the results in :bash:`workdir`:

    - :bash:`$workdir/obsvect/monitor.nc`: the output file with perturbed simulations

The example produces five days of simulations using a Gaussian model with arbitrary meteorological conditions,
fluxes distributed to show the text 'CIF', 30 stations randomly distributed over the domain.
You can play with these parameters to see the different possible outputs:

    - l. 104: :bash:`nstations`: you can reduce or increase the number of virtual sites
    - l. 203: :bash:`flx_txt`: you can change the text you want to be used for generating the fluxes; use short string to avoid having a text too pixelized


