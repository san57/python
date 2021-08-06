###################################
Standard pyCIF observation operator
###################################


.. role:: bash(code)
   :language: bash


Observation operators compute the operation
:math:`\mathbf{x} \rightarrow \mathcal{H}(\mathbf{x})` and its adjoint
if necessary. This operation relies on simulations by the chosen
numerical model, which are automatically initialized, articulated and
chained with each other by the observation operator.

Only one observation operator is coded at the moment, with the following
options:

-  :bash:`onlyinit` (optional): if True, the observation operator will only
   initialize model inputs; no simulation is run in this mode, which
   might lead to bugs if the observation operator is called by a
   function expecting some outputs; the recommended mode for this option
   is the :bash:`forward` execution mode

-  :bash:`autorestart` (optional): if True, the operator will check whether
   the numerical simulations were already computed; if it is the case,
   jumps to the last valid simulations; this is to be use to run a
   simulation that crashed; remove the last sub-directory to be
   re-computed; default is False

Therefore, the following Yaml paragraph is optional to initialize the
observation operator in default mode, but is required for specifying
additional options:

.. code-block:: Yaml

    obsoperator:
      plugin:
        name: standard
        version: std
      onlyinit: True
      autorestart: False

