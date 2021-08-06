#######################
Gaussian Toy Model
#######################


.. role:: bash(code)
   :language: bash

Description
-----------

These are the time series of meteorological parameters (wind speed, wind direction, stability class) needed to run the Gaussian Toy Model.


Configuration
-------------

The following Yaml paragraph is necessary for the configuration of the Gaussian Toy Model:

.. code-block:: Yaml

    meteo :
      plugin :
        name : dummy
        version : csv
      dirmeteo : /home/users/aberchet/CIF/data/dummy_gauss/
      filemeteo : meteo2.csv
      resolution: '1H'

The Yaml arguments are:

-  :bash:`dirmeteo`:
    directory where to find meteorological file

-  :bash:`filemeteo`:
    the meteo file in the directory

-  :bash:`resolution`:
    temporal resolution to be used to create meteo time series if :bash:`filemeteo` does not exist


