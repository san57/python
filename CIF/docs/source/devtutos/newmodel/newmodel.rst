############################
How to implement a new model
############################

.. role:: bash(code)
   :language: bash

In the CIF, a given numerical model has two main components:

    1. a black-box:
        it includes physical, chemical, transport processes, etc.;
        in most cases it is an executable compiled from Fortran or C sources;
        it computes outputs (concentration fields, footprints, etc.) from prescribed inputs
        (fluxes, meteorological fields, boundary conditions, etc.)
    2. a Python interface to the model:
        it produces the model inputs, reads the outputs and runs the executable properly.

The Python interface to a model is implemented as a pyCIF :bash:`model` class.

Prior to starting to implement your model, please make sure you have the following elements:

    1. a working version of your model (compilable sources), in forward mode, and if applicable in adjoint and tangent-linear modes

    2. a test directory with compatible inputs and outputs for your model

    3. an observation file you use to compare with your model


In the following, you will learn how to implement a model to pyCIF step-by-step.


.. toctree::
    :maxdepth: 3

    sources.rst
    register.rst
    monitor.rst
    obsvect.rst
    out2native.rst
