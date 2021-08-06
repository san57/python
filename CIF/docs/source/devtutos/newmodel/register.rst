#####################################
Create and register the plugin module
#####################################

Here, one should add the model to the list of model plugins recognized by the CIF.
The list of available plugins can be accessed :doc:`here<../../documentation/plugins/available>`

1. Go to the directory "pycif/plugins/models/"
2. Make an empty directory named according to your model.
3. Create an empty file called :bash:`__init__.py` to turn the folder into a Python module:

.. code-block:: bash

    cd CIF/pycif/plugins/models
    mkdir your/model/
    cd your/model
    touch __init__.py

4. Edit the register.py file and add the following lines:

.. code-block:: python

    import your-new-model
    Model.register_plugin('Your-Model-Name', 'Some-version', your-new-model)

When registering, pay attention to not replicate an already existing model!

5. Check that your model is available in pyCIF.
Your model should appear in the list when typing the following commands in python:

.. code-block:: python

    from pycif.utils.classes.baseclass import Plugin
    Plugin.print_registered()

returns:

.. code-block:: python

    List of all available plugins for each class:


     domain
    ---------------------------
        - CHIMERE, std
        - dummy, std
        - LMDZ, std


     chemistry
    ---------------------------
        - CHIMERE, gasJtab


     obsoperator
    ---------------------------
        - standard, std


     obsparser
    ---------------------------
        - WDCGG, std


     simulator
    ---------------------------
        - gausscost, std
        - dummy_txt, std


     controlvect
    ---------------------------
        - standard, std


     measurements
    ---------------------------
        - random, std
        - standard, std


     meteo
    ---------------------------
        - dummy, csv
        - CHIMERE, std
        - LMDZ, mass-fluxes


     platform
    ---------------------------
        - LSCE, obelix


     obsvect
    ---------------------------
        - standard, std


     mode
    ---------------------------
        - analytic, std
        - adj-tl_test, std
        - 4dvar, std
        - footprint, std
        - forward, std
        - post-proc, std


     model
    ---------------------------
        - dummy, std
        - LMDZ, std
        - CHIMERE, std


     fluxes
    ---------------------------
        - LMDZ, sflx
        - dummy, txt
        - dummy, nc
        - LMDZ, bin
        - CHIMERE, AEMISSIONS


     minimizer
    ---------------------------
        - congrad, std
        - M1QN3, std
