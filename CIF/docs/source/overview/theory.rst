Theoretical framework
---------------------

.. role:: raw-math(raw)
    :format: latex html

The CIF should be compatible with analytical inversions, variational inversions, as well as Ensemble Kalman Filters (EnKFs) as the main methods used in the community.
All these methods rely on the classical Bayesian inversion framework with Gaussian assumptions.
The Gaussian Bayesian formulation of the inversion problem consists in computing the following probability density function (pdf):

.. math::
    :label: bayes

    \begin{equation}
    p(\mathbf{x} | \mathbf{y}^\textrm{o}, \mathbf{x}^\textrm{b}) \sim \mathcal{N}(\mathbf{x}^\textrm{a}, \mathbf{P}^\textrm{a})
    \end{equation}


with :math:`\mathbf{y}^\textrm{o}` the observation network,
:math:`\mathbf{x}^\textrm{b}` rendering the prior knowledge on variables
to optimize (most of the time fluxes in our case, but also concentration
fields in some configurations).

Analytical inversions
^^^^^^^^^^^^^^^^^^^^^

When the observation operator is linear, :math:`\mathcal{H}` can be
fully described by its Jacobian matrix :math:`\mathbf{H}`, and
conversely its adjoint :math:`\mathcal{H}^*` by the transpose of the
Jacobian :math:`\mathbf{H}^\textrm{T}`. Thus,
:math:`\mathbf{x}^\textrm{a}` and :math:`\mathbf{P}^\textrm{a}` can be
explicitly written as:

.. math::
    :label: analytical

    \begin{equation}
    \left\{
    \begin{array}{rclcl}
    \mathbf{x}^\textrm{a} & = & \mathbf{x}^\textrm{b} + \mathbf{K}(\mathbf{y}^\textrm{o} - \mathbf{H}\mathbf{x}^\textrm{b}) \\
    \mathbf{P}^\textrm{a} & = & \mathbf{P}^\textrm{b} - \mathbf{K}\mathbf{H}\mathbf{P}^\textrm{b}
    \end{array}
    \right.
    \textrm{ or }
    \left\{
    \begin{array}{rclcl}
    \mathbf{x}^\textrm{a} & = & \mathbf{x}^\textrm{b} + \left[\mathbf{H}^\textrm{T}\mathbf{R}^{-1}\mathbf{H} + (\mathbf{P}^\textrm{b})^{-1}\right]^{-1} \mathbf{H}^\textrm{T} (\mathbf{y}^\textrm{o} - \mathbf{H}\mathbf{x}^\textrm{b}) \\
    \mathbf{P}^\textrm{a} & = & \left[\mathbf{H}^\textrm{T}\mathbf{R}^{-1}\mathbf{H} + (\mathbf{P}^\textrm{b})^{-1}\right]^{-1} \\
    \end{array}
    \right.
    \end{equation}

with :math:`\mathbf{K}` the Kalman gain matrix:
:math:`\mathbf{K} = \mathbf{P}^\textrm{b}\mathbf{H}^\textrm{T}(\mathbf{R}+\mathbf{H}\mathbf{P}^\textrm{b}\mathbf{H}^\textrm{T})^{-1}`


The formulation with the Kalman gain matrix is limited by the inversion
of a matrix of dimension dim(\ :math:`\mathcal{Y}`), the observation
space dimension, while the other formulation is limited by the dimension
of the control space. Due to the very high dimensions for both the
observation and the control spaces in most inversion applications, the
explicit computation of Eq. :eq:`analytical` with
matrix products and inverses is not computationally feasible. For this
reason, smart adaptations on the inversion framework (including
approximations and numerical solvers) are necessary to tackle the
problem.

Variational inversions
^^^^^^^^^^^^^^^^^^^^^^

One possible way to avoid the dimension issue is the variational
approach. Computing the normal distribution in Eq. :eq:`bayes` is equivalent to finding the
minimum of the cost function:

.. math::
    :label: variational

    \begin{equation}
    J(\mathbf{x}) =
    \frac{1}{2} (\mathbf{x} - \mathbf{x}^\textrm{b})^\textrm{T}
    (\mathbf{P}^\textrm{b})^{-1} (\mathbf{x} - \mathbf{x}^\textrm{b})
    + \frac{1}{2} (\mathcal{H}(\mathbf{x}) - \mathbf{y}^\textrm{o})^\textrm{T}
    \mathbf{R}^{-1}(\mathcal{H}(\mathbf{x}) - \mathbf{y}^\textrm{o})
    \end{equation}

In variational inversions, the minimum of the cost function in Eq. :eq:`variational` is numerically computed using
a quasi-Newtonian descending algorithm based on the gradient of the cost
function:

.. math::
    :label: gradient

    \begin{equation}
    \nabla J_\mathbf{x} = (\mathbf{P}^\textrm{b})^{-1} (\mathbf{x} - \mathbf{x}^\textrm{b})
    + \mathcal{H}^*\left[\mathbf{R}^{-1}(\mathcal{H}(\mathbf{x}) - \mathbf{y}^\textrm{o})\right]
    \end{equation}

Quasi-Newtonian methods are a group of algorithms designed to compute
the minimum of a function. In the community, one example of
quasi-Newtonian algorithms commonly used is M1QN3 (`Gilbert and
Lemaréchal,
1989 <https://link.springer.com/article/10.1007/BF01589113>`__). In
general quasi-Newtonian methods require an initial regularization of
:math:`\mathbf{x}`, the vector to be optimized, for better effileiency.
In atmospheric inversion, such a regularization is generally made by
optimizing
:math:`\mathbf{\chi} = (\mathbf{P}^\textrm{b})^{-1/2} (\mathbf{x} - \mathbf{x}^\textrm{b})`
instead of
:math:`\mathbf{x}`. Although more optimal regularizations can be chosen,
the minimization of the equations with :math:`\mathbf{\chi}` is
preferred for its simplifying the equation to solve. This transformation
translates in Eq. :eq:`gradient` as follows:

.. math::
    :label: gradientchi

    \begin{equation}
    \nabla J_\mathbf{\chi} = \chi
    + (\mathbf{P}^\textrm{b})^{1/2}\mathcal{H}^*\left[\mathbf{R}^{-1}(\mathcal{H}(\mathbf{x}) - \mathbf{y}^\textrm{o})\right]
    \end{equation}

Ensemble Kalman Filters
^^^^^^^^^^^^^^^^^^^^^^^

In EnKFs, such as presented in e.g., `Peters et al. (2005) <https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/2005JD006157>`__,
the issue of the high dimension in the system of Equations :eq:`analytical` is avoided using two main procedures:

- observations are assimilated sequentially in the system to reduce the dimension of the observation space, making it possible to compute matrix products and inverses
- covariance matrices are approximated with a Monte Carlo ensemble of possible control vectors:

.. math::
    :label: enkf

    \begin{equation}
    \label{eq:enkf}
    \left\{
    \begin{array}{rcl}
    \mathbf{H}\mathbf{P}^\textrm{b}\mathbf{H}^\textrm{T} & \simeq & \frac{1}{N-1}(\mathcal{H}(\mathbf{x}_1), \mathcal{H}(\mathbf{x}_2), ..., \mathcal{H}(\mathbf{x}_N))\cdot(\mathcal{H}(\mathbf{x}_1), \mathcal{H}(\mathbf{x}_2), ..., \mathcal{H}(\mathbf{x}_N))^\textrm{T} \\
    \mathbf{P}^\textrm{b}\mathbf{H}^\textrm{T} & \simeq & \frac{1}{N-1}(\mathbf{x}_1, \mathbf{x}_2, ..., \mathbf{x}_N)\cdot(\mathcal{H}(\mathbf{x}_1), \mathcal{H}(\mathbf{x}_2), ..., \mathcal{H}(\mathbf{x}_N))^\textrm{T} \\
    \end{array}
    \right.
    \end{equation}


