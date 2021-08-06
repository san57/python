File name format
----------------

.. role:: bash(code)
   :language: bash


Various types of files are necessary to run pyCIF. pyCIF accepts generic
formats when file names are dependent on, e.g., the date of the
simulation. In this case, the generic date-dependent file name should be
given as a string compatible with the Python datetime function,
:bash:`strptime` as detailed
`here <https://docs.python.org/2/library/datetime.html#strftime-and-strptime-behavior>`__.

For instance, if your meteo file names are function of the current date,
the following syntax will be required in the Yaml file:

.. code-block:: yaml

    mymeteo: /mymeteo_folder/mymeteo_%Y%m%d_%H%M.nc

pyCIF will automatically recognize a date format as YYYYMMDD_HHMM, i.e.
20140327_1212 corresponds to March 27 :sup:`th`, 2014
at 12:12.
