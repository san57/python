
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
    # How to build your observation vector and observation uncertainties if needed
    # Also projects information from the observation to the model space
    # - file_obsvect: observation vector from previous simulations
    obsvect:
      plugin:
        name: standard
        version: std
      file_obsvect : /your/monitor/file
      dump: True
    #####################################################################
    #####################################################################


