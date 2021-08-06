###################
CONGRAD
###################


.. role:: bash(code)
   :language: bash

Description
-----------

CONGRAD follows a conjugate gradient method combined with a Lanczos
algorithm. It can be used to solve positive definite quadratic functions
or to solve linear systems involving a symmetric positive definite
matrix. It was originally developed at ECMWF (`Fisher,
1998 <https///www.ecmwf.int/sites/default/files/elibrary/1998/9400-minimization-algorithms-variational-data-assimilation.pdf>`__)
and adapted in Python in 2004.

The algorithm requires the following arguments:

-  :bash:`maxiter`: maximum number of iterations; mandatory

-  :bash:`zreduc`: required reduction in gradient norm; optional; default is
   1e-15

-  :bash:`pevbnd`: Accuracy required of approximate eigenvectors; optional;
   default is 0.01

-  :bash:`kverbose`: verbose level, from 0 to 2; optional; default is 1

-  :bash:`ldsolve`: minimize the cost function; optional; default is 1; for
   expert users only

A Yaml template presents as follows:

.. code-block:: Yaml

    :::Yaml
    minimizer :
      plugin:
        name: congrad
        version: std
      simulator:
        plugin:
          name: gausscost
          version: std
      kverbose: 1

Requirements
------------

CONGRAD requires the following plugins to be executed properly: 

1. a :doc:`simulator </configuration/simulators/index>` to compute function
to minimize; optional: default is (:bash:`gausscost`, :bash:`std`)
