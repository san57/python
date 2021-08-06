######################
Variational inversions
######################

.. role:: bash(code)
   :language: bash


Variational inversion
---------------------

Description
-----------

This mode computes a full variational inversion based on measurements, a
corresponding observation vector, and a control vector. This mode needs
a minimizer to be computed, which itself requires a so-called
''simulator'', i.e. the function to minimize (at the moment, the default
function is the Bayesian Gaussian cost function).

Details on the minimizers are given :doc:`here </configuration/minimizers/index>`.

The Yaml paragraph corresponding to a variational run is as follows:

.. code-block:: Yaml

    mode:
      plugin:
        name: 4dvar
        version: std
      minimizer :
        plugin:
          name: congrad
          version: std
        simulator:
          plugin:
            name: gausscost
            version: std
        maxiter: 10
        epsg: 0.03
        df1: 0.01
        kverbose: 2

Requirements
------------

A variational inversion requires the following plugins to be executed
properly:

1. a :doc:`control vector </configuration/controlvect/index>`
to define the control vector shape and corresponding operations;
mandatory

2. an :doc:`observation vector </configuration/obsvect/index>` to define observations to
be compared with; mandatory

3. an :doc:`minimizer </configuration/minimizers/index>` to optimize the
function; optional; default is (:bash:`m1qn3`, :bash:`std`)

The following plugins are indirectly needed to compute a variational
inversion:

1. a :doc:`numerical model </configuration/models/index>`; mandatory

2. an :doc:`observation operator </configuration/obsoper/index>` to compute
:math:`\mathbf{x} \rightarrow \mathcal{H}(\mathbf{x})` and its adjoint;
optional: default is (:bash:`standard`, :bash:`std`)

Example
-------




