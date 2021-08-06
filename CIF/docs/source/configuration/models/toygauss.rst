#######################
Gaussian Toy Model
#######################


.. role:: bash(code)
   :language: bash

Description
-----------

The Gaussian Toy Model was designed to quickly test pycif.
It is a simple Gaussian model explicitly computing the source-receptor relationship based on the equations of the Gaussian plume.
Meteorological fields (wind fields and stability class) are assumed constant on the entire domain.
The parameters in the Gaussian plume equations are dependent of the wind speed and the stability class

Configuration
-------------

The following Yaml paragraph is necessary for the configuration of the Gaussian Toy Model:

.. code-block:: Yaml

    model :
        plugin:
        name    : dummy
        version : std

        # H matrix
        file_pg : /home/users/aberchet/CIF/model_sources/dummy_gauss/Pasquill-Gifford.txt

The Yaml arguments are:

-  :bash:`file_pg`:
    the reference Pasquill-Gifford definition of stability class and corresponding wind speed ranges;
    this file is provided with the CIF in the :bash:`data` folder

Requirements
------------









