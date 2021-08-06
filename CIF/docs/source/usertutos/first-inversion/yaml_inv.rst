.. code-block:: yaml
    :linenos:

    #####################
    # PYCIF config file #
    #####################

    # Define here all parameters for PYCIF following Yaml syntax
    # For details on Yaml syntax, please see:
    # http://docs.ansible.com/ansible/latest/YAMLSyntax.html

    # For non-specified parameters, one can either comment the corresponding line
    # or give the value "null"

    # Some parameters are switches that need to be activated or not.
    # PYCIF expects a boolean True/False
    # Yaml accepts the following syntax to be converted to booleans:
    # y|Y|yes|Yes|YES|n|N|no|No|NO|true|True|TRUE|false|False|FALSE|on|On|ON|off|Off|OFF
    # This offers flexibility, but one should stick to a style for consistency

    #####################################################################
    #####################################################################
    # PYCIF parameters

    # Verbose level
    # 1 = Basic information
    # 2 = Debugging
    verbose : 1

    # Log file (to be saved in $wordkir)
    logfile: pycif.logtest

    # Execution directory
    workdir : ${HOME}/CIF_dummy_adjtl_test/

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
    # To-Do for initializing PYCIF
    # Can be commented if running PYCIF directly
    #todo_init: [measurements]
    #####################################################################
    #####################################################################


    #####################################################################
    #####################################################################
    # Running mode for PYCIF
    mode:
      plugin:
        name: 4dvar
        version: std
      minimizer:
        plugin:
          name: congrad
          version: std
        simulator:
          plugin:
            name: gausscost
            version: std
          reload_from_previous: True
        maxiter: 20
        epsg: 0.02
        df1: 0.01
    #####################################################################
    #####################################################################


    #####################################################################
    #####################################################################
    # Observation operator
    # Used to run the model and translates information from model/measurement
    # spaces to control/observation spaces
    obsoperator:
      plugin:
        name: standard
        version: std
      autorestart: True
    #####################################################################
    #####################################################################


    #####################################################################
    #####################################################################
    # Transport model
    # Accepts any model registered in pycif.models
    model :
      plugin:
        name    : dummy
        version : std

      # H matrix
      file_pg : ~/cif/model_sources/dummy_gauss/Pasquill-Gifford.txt

      # Chemical scheme
      chemistry:
        acspecies:
          CH4:
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
      file_monitor :  ~/CIF_dummy_adjtl_test/monitor_reference.nc
      dump_type : nc
      species :
        # List of tracers to include in the monitor file
        # For each tracer observations, please specify:
        # - provider = list of observation providers
        # - format = list of format types
        # - dir_obs = directory where to find observation files
        # - err_obs (optional) = uncertainty on observations
        # - depos (optional) = surface deposition file
        # - rescale (optional) = true to change observation scale if necessary
        # - na_values (optional) = invalid values to exclude. Default is -999
        # - default_unit (optional) = basic unit for reporting the tracer
        #                             everything is then converted to ppm
        # - dump (optional) = dump to a monitor file. Default is True
        # For tracer fluxes, please specify:
        # - dir_flx: directory where reference fluxes are already computed
        # - file_flx: file format to recover fluxes from
        CH4 :
          frequency: '1H'
          nstations: 50
          duration: '2H22min'
          random_subperiod_shift: True
          zmax: 200
        MCF :
          frequency: '1H'
          nstations: 20
          duration: '5H22min'
          random_subperiod_shift: True
          zmax: 200

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
      file_obsvect : ~/CIF_dummy_adjtl_test/monitor_reference.nc
      dump_type: nc
      transform:
        timeavg:
          plugin:
            name: timeavg
            version: std
            type: transform
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
    statevect:
      plugin:
        name: standard
        version: std
      components:
        fluxes:
          parameters:
            CH4 :
              plugin:
                name: 'dummy'
                version: 'txt'
                type: 'fluxes'
              hresol : hpixels
              type : physical
              errtype : avg
              xb_value : 1
              err : 1
              period : '36H'
              flx_text: 'CIF'
              dir: ~/CIF_dummy_adjtl_test/statevect/
              file: flx_real.txt
              hcorrelations :
                landsea: False
                dump_hcorr : True
                dircorrel : ~/CIF_dummy_adjtl_test/
                sigma: 1000
              tcorrelations :
                sigma_t: 200
        meteo :
          plugin :
            name : dummy
            version : csv
            type: meteo
          dir: ~/CIF/data/dummy_gauss/
          file : meteo2.csv
          resolution: '1H'
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
      nlon: 60
      ymin: 0
      ymax: 20000
      nlat: 30
    #####################################################################
    #####################################################################

