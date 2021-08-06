#####################
Toy Gaussian model
#####################

This toy model is implemented in the CIF as a testing utility.
It can easily be run on office laptop with low computational requirements.
Thus, features in pycif can be explored quickly.

It is also proposed as a minimal example of a fully implemented model in pycif.
We recommend developers to explore its code as a template for new model implementation in pycif.

***********************
Gaussian plume equation
***********************

The Toy Gaussian Model follows the standard version of the Gaussian plume equation
based on Pasquill-Gifford stability classification:


.. math::
    :label: gaussplume

    \begin{equation}
    \frac{1}{2 \pi \sigma_y \sigma_z \bar{u}}
    \exp\left(-\frac{y^2}{\sigma_y ^2}\right)
    \exp\left(-\frac{z^2}{\sigma_z ^2}\right)
    \end{equation}

with


.. math::
    :label: gaussplume

    \begin{equation}
    \left\{
    \begin{array}{rcl}
    \sigma_z &=& ax^b   \\
    \sigma_y &=& | 465.11628 x \tan(0.017653293 (c - d \ln x)) |
    \end{array}
    \right.
    \end{equation}

x is the distance between the source and receptor points, (a, b, c, d) are parameters depending on the stability class.



