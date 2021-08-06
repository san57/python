######################
M1QN3
######################


.. role:: bash(code)
   :language: bash

Description
-----------

M1QN3 is a quai-Newtonian optimization algorithm originally developed at
`INRIA <https///who.rocq.inria.fr/Jean-Charles.Gilbert/modulopt/optimization-routines/m1qn3/m1qn3.html>`__
to minimize functions with a very high number of variables. The
algorithm is described in `Gilbert and Lemaréchal,
1989 <https///link.springer.com/article/10.1007/BF01589113>`__. The
original code was written in Fortran; F. Chevallier translated it to
Python in 2005.

The algorithm requires the following arguments:

- :bash:`niter`: maximum number of iterations

- :bash:`nsim`: maximum number of simulations

- :bash:`maxiter`: maximum number of iterations; if one of the two previous
   is missing,:bash:`niter` is set to:bash:`maxiter` and:bash:`nsim` to 2 times
  :bash:`maxiter`

- :bash:`dxmin`: absolute precision on x; optional; default is 1e-20

- :bash:`df1`: expected decrease for f; optional; default is 0.01

- :bash:`epsg`: relative precision on the gradient; optional; default is
   1e-20

- :bash:`mode`: mode more M1QN3; optional; default is 0; for expert users
   only

- :bash:`m`: number of updates; optional; default is 5; for expert users
   only

A Yaml template presents as follows:

.. code-block:: Yaml

    minimizer :
      plugin:
        name: M1QN3
        version: std
      simulator:
        plugin:
          name: gausscost
          version: std
      maxiter: 10
      epsg: 0.03
      df1: 0.01

Requirements
------------

M1QN3 requires the following plugins to be executed properly: 1. a
:doc:`simulator </configuration/simulators/index>` to compute function
to minimize; optional: default is (:bash:`gausscost`,:bash:`std`)
