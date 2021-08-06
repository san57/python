#####################################
How to run a first inversion
#####################################

.. role:: bash(code)
   :language: bash

After running a :doc:`first forward run<../first-simu/index>` with perturbed outputs from known emissions with the :doc:`the Toy Gaussian Model</documentation/plugins/models/toy-gaussian>`,
it is possible to carry out an academic inversion using these "true" observations.
Below are the steps to run a variational inversion with the minimizer :bash:`congrad`:

1. Set up a Yaml configuration:

    - copy the content of the reference Yaml file available :doc:`here<yaml_inv>` to a directory of your choice
    - modify paths to fit your own configuration:
        - l. 34: :bash:`workdir`: where outputs will be generated
        - l. 67: :bash:`file_pg`: Pasquil-Gifford data, provided with CIF codes
        - l. 122: :bash:`file_obsvect`: observations fitted to the model grid; should point to the :bash:`monitor.nc` file produced in our :doc:`first forward run<../first-simu/index>`
        - l. 152: :bash:`dir_meteo`: directory to example meteo files, provided with CIF codes
        - l. 194: :bash:`dir`: directory where fluxes will be generated
        - l. 217: :bash:`flx_text`: the text used to generate prior (should be different than the text used in the reference forward run)

2. Run pyCIF based on your Yaml file:

    .. code-block:: bash

        python -m pycif path-to-your-yaml

3. Explore the results in :bash:`workdir`:

    - :bash:`$workdir/controlvect/fluxes/controlvect_fluxes_CH4.nc`: the optimized fluxes
    - :bash:`$workdir/controlvect/simulator/cost.txt`: the value of the cost function at each iteration
    - :bash:`$workdir/controlvect/simulator/gradcost.txt`: the norm of the gradient of the cost function at each iteration

You just carried out an inversion using the minimizer :bash:`congrad`.

.. note::
    It is very likely that your inversion did not succeed in capturing the target text to retrieve.
    You can play with the number of stations to test how many are needed to retrieve the 'truth' with a very bad prior.

