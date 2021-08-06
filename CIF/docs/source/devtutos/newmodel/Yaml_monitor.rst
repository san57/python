
.. code-block:: yaml

    #####################
    # pyCIF config file #
    #####################

    # Define here all parameters for pyCIF following Yaml syntax
    # For details on Yaml syntax, please see:
    # http://docs.ansible.com/ansible/latest/YAMLSyntax.html

    # For non-specified parameters, one can either comment the corresponding line
    # or give the value "null"

    # Some parameters are switches that need to be activated or not.
    # pyCIF expects a boolean True/False
    # Yaml accepts the following syntax to be converted to booleans:
    # y|Y|yes|Yes|YES|n|N|no|No|NO|true|True|TRUE|false|False|FALSE|on|On|ON|off|Off|OFF
    # This offers flexibility, but one should stick to a style for consistency

    #####################################################################
    #####################################################################
    # pyCIF parameters

    # Verbose level
    # 1 = Basic information
    # 2 = Debugging
    verbose : 1

    # Log file (to be saved in $wordkir)
    logfile: pyCIF.logtest

    # Execution directory
    workdir : /some_working_directory/

    # Initial date
    # Use the following compatible format for the date:
    # YYYY-mm-dd or YYYY-mm-dd HH:mm
    datei : 2010-01-01

    # End date (same as for the initial date)
    datef : 2010-01-05
    #####################################################################
    #####################################################################


    #####################################################################
    #####################################################################
    # Measurements to account for
    # Main keys are:
    # - infos: infos on tracers to include in the simulation
    # - file_monitor: standard txt file including all observations
    measurements :
      plugin:
        name: random
        version: std

     # File where to save data, if does not exist. Reads from there if exists
      file_monitor :  /some_directory/monitor_reference.nc
      dump_type : nc
      species :
        # List of tracers to include in the monitor file
        # For each tracer observations, please specify:
        #   - nstations = number of stations to randomly sample
        #   - duration = duration of each individual data point
        #   - frequency = for each station, one obs every XX frequency
        #   - random_subperiod_shift = randomly shifts observations within frequency intervals
        #   - zmax = maximal vertical altitude to use for stations
        some_species :
          frequency: '1H'
          nstations: 30
          duration: '1H'
          random_subperiod_shift: True
          zmax: 100
    #####################################################################
    #####################################################################


    #####################################################################
    #####################################################################
    # Domain definition
    domain :
      plugin :
        name : dummy
        version : std
      xmin: -180
      xmax: 180
      nlon: 360
      ymin: -90
      ymax: 90
      nlat: 180
    #####################################################################
    #####################################################################
