#################################
Standard pyCIF observation vector
#################################


.. role:: bash(code)
   :language: bash


Observation vectors store observation points (not necessary exactly raw
measurements) and include related operations, such as matrix product
with observation uncertainties.

Only a standard observation vector is implemented at the moment.

Configuration
-------------

The observation vector loads functions, but also observation data. This
is done with the following arguments:

- :bash:`file_obsvect`: the file with the observation data (see
    :doc:`here </configuration/file_formats>` for
    further details on the format); optional; if not specified, tries
    computing measurements as set-up :doc:`here </configuration/measurements/index>`; if computation
    fails (due to, e.g., empty measurements), passes an empty data store
    for later operations

- :bash:`dump`: dumps newly computed observation vector to:bash:`file_obsvect`
   if True; optional; default is True

- :bash:`dump_type` (optional) = dump monitor either as a NetCDF (:bash:`nc`)
   or a CSV file (:bash:`csv`); default is :bash:`csv`

- :bash:`file_statlev`: file specifying station vertical in the level
   (should soon be generalized)

- :bash:`wfunc`: weighting functions for satellites; optional; default is
   False


Requirements
------------

The observation vector requires the following plugins to be executed
properly: 

1. :doc:`measurements </configuration/measurements/index>`
storing un-processed observations; optional: default is an empty
datastore

2. a :doc:`model </configuration/models/index>` to define the
control vector shape and corresponding operations; mandatory
