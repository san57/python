#############################################
Test of the adjoint and of the linear-tangent
#############################################

.. role:: bash(code)
   :language: bash


Description
-----------

For every inversion case study, it is necessary to compute a test of the
adjoint before running a full variational inversion. It allows to detect
any error specifi to your test case and that were never met before. If
the computation costs of the inversion set-up are very high, it is
possible to execute a adjoint test with a simplified set-up, but it is
important to include configurations that can generate errors (month of
February, several months/years in a row, etc.).

This test consists in checking that:
:math:`\langle \mathcal{H}(\delta \mathbf{x}) \vert \mathcal{H}(\delta \mathbf{x}) \rangle = \langle \delta \mathbf{x} \vert \mathcal{H}^* \circ \mathcal{H}(\delta \mathbf{x}) \rangle`.

:math:`\mathcal{H}(\delta \mathbf{x})` is computed with the tangent
linear model, and then
:math:`\mathcal{H}^*\circ \mathcal{H}(\delta \mathbf{x})` is calculated
with the adjoint.

In principle, the two terms of the equation are never exactly equal.
Nevertheless, the difference should never exceed several times the
machine accuracy.

Configuration
-------------

The following Yaml template is valid for the adjoint test:

.. code-block:: Yaml

    mode:
      plugin:
        name: adj-tl_test
        version: std

      # Magnitude of the increment for dx
      increments: 0.2

      # (optional) Type of increments: 
      # - 'cst' for all equal to 'increments'; 
      # - 'rand' for a centered normal distribution with magnitude 'increments'
      # Default is 'cst' if not specified
      incrmode: cst

      # Testing space. Should be one of: control, chi
      testspace: chi
      

The Yaml arguments are:

-  :bash:`increments`: determines the value of :math:`\delta \mathbf{x}`,
   which should be taken within a reasonable range around 0

-  :bash:`incrmode`: the type of increments

   -  should be one of: :bash:`cst`, :bash:`rand`
   -  the default value is :bash:`cst` if not specified
   -  if :bash:`cst` all elements of :math:`\delta \mathbf{x}` are set to
      the same value
   -  if :bash:`rand` elements of :math:`\delta \mathbf{x}` follow a
      centered normal distribution with standard deviation
      :bash:`increments`
   -  it is recommended to execute an adjoint test with both mode to
      make sure that some errors do not compensate each other

-  :bash:`testspace`: testing space

   -  should be one of: :bash:`control`, :bash:`chi`; default is :bash:`control`
      when not specified
   -  if :bash:`control` the test is applied to :math:`\delta \mathbf{x}` in
      the control space
   -  if :bash:`chi` the test is computed in the minimization space on
      :math:`\delta \mathbf{\chi}`; it allows validating operations at
      the control space level, e.g.,
      :math:`\delta \mathbf{x} \rightarrow \mathbf{B}^{1/2} \delta \mathbf{x}`

Requirements
------------

The adjoint test requires the following plugins to be executed properly:

1. an :doc:`observation operator </configuration/obsoper/index>` to
compute :math:`\mathbf{x} \rightarrow \mathcal{H}(\mathbf{x})` and its
adjoint; optional: default is (`standard`, :bash:`std`) 

2. a :doc:`control vector </configuration/controlvect/index>` to define the control
vector shape and corresponding operations; mandatory 

3. an :doc:`observation vector </configuration/obsvect/index>` to define observations to
be compared with; mandatory

The following plugins are indirectly needed to compute the test of the adjoint,
through the observation operator (see details
:doc:`here </configuration/obsoper/index>`):

1. a :doc:`numerical model </configuration/models/index>`; mandatory


Example
-------
