
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
    workdir : ~/CIF_dummy_test/

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
    # Running mode for pyCIF
    mode:
      plugin:
        name: forward
        version: std
      perturb_obsvect: True
      obserror: 0.01
    #####################################################################
    #####################################################################


    #####################################################################
    #####################################################################
    # Transport model
    # Accepts any model registered in pycif.plugins.models
    model :
      plugin:
        name    : dummy
        version : std

      # H matrix
      file_pg : ~/CIF/model_sources/dummy_gauss/Pasquill-Gifford.txt
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
      file_monitor :  ~/CIF_dummy_test/monitor_reference.nc
      dump_type : nc
      species :
        CH4 :
          frequency: '1H'
          nstations: 30
          duration: '1H'
          random_subperiod_shift: True
          zmax: 100

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
      file_obsvect : ~/CIF_dummy_test/monitor_reference.nc
      dump_type: nc
    #####################################################################
    #####################################################################


    #####################################################################
    #####################################################################
    # Arguments to define the state vector
    # These are specifi to the state vector and inversion method you chose
    # Please refer to the documentation to know what to include here
    # For the standard LMDZ, include the following:
    # - filelsm: land-sea mask (must be consistent with LMDZ grid)
    # - correl: Use correlation or not in B
    # - dircorrel: path to pre-computed correlations
    # - sigma_land: spatial correlation length for prior errors over land (km)
    # - sigma_sea: spatial correlation length for prior errors over ocean (km)
    # - tracers: list of tracers to put in the state vector (with definition arguments):
    #     - calcstd: calculate global standard deviation
    #     - hresol: resolution at which fields are scaled
    #            (choice = bands,regions,pixels;
    #             if regions, provide a netcdf file fileregion
    #             if bands, define a list of latitudes as band limits (n+1 for n bands)
    #     - periodflux: period of variation for increments within a month (days)
    #     - glob_err (optional) = uncertainty on global budget
    controlvect:
      plugin:
        name: standard
        version: std
      save_out_netcdf: True
      components:
        fluxes:
          parameters:
            CH4 :
              plugin:
                name: 'dummy'
                version: 'txt'
              hresol : hpixels
              type : physical
              errtype : max
              err : 1
              period : '1D'
              dir: ~/CIF_dummy_test/controlvect/
              file: flx_real.txt
              hcorrelations :
                landsea: False
                dump_hcorr : True
                dircorrel : ~/CIF_dummy_test/controlvect/
                sigma: 1000
              tcorrelations :
                sigma_t: 5
              flx_text: 'CIF'
    #####################################################################
    #####################################################################


    #####################################################################
    #####################################################################
    # Domain definition
    domain :
      plugin :
        name : dummy
        version : std
      xmin: 0
      xmax: 25000
      nlon: 30
      ymin: 0
      ymax: 20000
      nlat: 15
    #####################################################################
    #####################################################################


    #####################################################################
    #####################################################################
    # Meteo definition
    meteo :
      plugin :
        name : dummy
        version : csv
      dirmeteo : ~/CIF/data/dummy_gauss/
      filemeteo : meteo.csv
    #####################################################################
    #####################################################################

