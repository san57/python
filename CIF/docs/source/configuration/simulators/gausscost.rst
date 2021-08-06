######################
Gaussian cost function
######################



.. role:: bash(code)
   :language: bash


Description
-----------

The Bayesian cost function is the basic quadratic form derived from the
Gaussian assumption in the application of the Bayesian theorem for data
assimilation as described :doc:`here </overview/index>`.
The defining equation is:

.. math::

   J(\mathbf{x}) = \frac{1}{2} (\mathbf{x} - \mathbf{x}^\textrm{b})^\textrm{T} (\mathbf{P}^\textrm{b})^{-1} (\mathbf{x} - \mathbf{x}^\textrm{b}) 
   + \frac{1}{2} (\mathcal{H}(\mathbf{x}) - \mathbf{y}^\textrm{o})^\textrm{T}\mathbf{R}^{-1}(\mathcal{H}(\mathbf{x}) - \mathbf{y}^\textrm{o})

As a quadratic form on a finite dimension space, it necessary has a
global minimum that can be fund with available :doc:`minimizers </configuration/minimizers/index>` (which do not
necessarily guarantee that the computed minimum is the global minimum).

The basic Gaussian cost function does not require specifi configuration
parameters and the corresponding Yaml paragraph is simply: A Yaml
template presents as follows:

.. code-block:: Yaml

    simulator:
      plugin:
        name: gausscost
        version: std

Requirements
------------

The Bayesian Gaussian cost function requires the following plugins to be
executed properly: 

1. a :doc:`control vector </configuration/controlvect/index>` to define the control
vector shape and corresponding operations; mandatory 

2. an :doc:`observation vector </configuration/obsvect/index>` to define observations to
be compared with; mandatory 

3. an :doc:`observation operator </configuration/obsoper/index>` to compute
:math:`\mathbf{x} \rightarrow \mathcal{H}(\mathbf{x})` and its adjoint;
optional: default is (:bash:`standard`, :bash:`std`)

The following plugins are indirectly needed to compute a variational
inversion: 

1. a :doc:`numerical model </configuration/models/index>`; mandatory
