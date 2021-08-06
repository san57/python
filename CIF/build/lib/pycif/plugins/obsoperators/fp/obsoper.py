import pandas as pd
import datetime
import numpy as np
import os

from pycif.utils import path
from pycif.utils.check import verbose
from pycif.utils.datastores.dump import dump_datastore, read_datastore
from .check import check_inputs


def obsoper(self, inputs, mode,
            run_id=0, datei=datetime.datetime(1979, 1, 1),
            datef=datetime.datetime(2100, 1, 1),
            workdir='./',
            reload_results=False,
            **kwargs):
    """The observation operator.
    This function maps information from the control space to the observation
    space
    
    This version of 'obsoper' was  developed for use with FLEXPART
    backward runs.

    Gets concentrations/mixing ratios from FLEXPART flux sensitivities.

    For now assumes FLEXPART runs are stored in the structure used by
    FLEXINVERT+,
    .i.e, /station_id/YYYYMM/
    
    if FWD
    Turns back model-equivalents into observation space
    Hx(native) -> Hx(obs)

    Generates departures and multiply by R-1
    -> Hx - y0, R^(-1).(Hx - y0)

    Calculates gradient in observation space
       Jo'(p) = H^TR^-1(H(p) - y)

    
    Args:
        inputs: can be a control vector or an observation vector, depending
                on the mode
        mode (str): the running mode; always 'fwd' for the flexpart version
        run_id (int): the ID number of the current run; is used to define
                the current sub-directory
        datei, datef (datetime.datetime): start and end dates of the
                simulation window
        workdir (str): the parent directory of the simulation to run
        reload_results (bool): look for results from pre-computed simulation
                if any
    
    Returns:
        observation or control vector depending on the running mode
    """
    
    # Check that inputs are consistent with the mode to run
    check_inputs(inputs, mode)
    
    # If true, do not run simulations at this point
    read_only = getattr(self, 'read_only', True)
    print("read_only", read_only)
    
    # Initializing modules and variables from the setup
    model = self.model
    controlvect = self.controlvect
    obsvect = self.obsvect
    subsimu_dates = model.subsimu_dates
    
    if mode == 'fwd':
        controlvect = inputs
        obsvect.datastore.loc[:, 'sim'] = np.nan
    
    # Various initializations 
    
    # Get flexpart headers
    subdir = subsimu_dates[0].strftime("%Y%m")
    
    fp_header_glob = model.utils.flexpart_header.Flexpartheader()
    fp_header_glob.read_header(
        os.path.join(model.run_dir_glob,
                     obsvect.datastore.head(1)['station'][0].upper(),
                     subdir, 'header'))

    fp_header_nest = None
    if model.plugin.nested:
        fp_header_nest = model.utils.flexpart_header.Flexpartheader()
        fp_header_nest.read_header(
            os.path.join(model.run_dir_nest,
                         obsvect.datastore.head(1)['station'][0].upper(),
                         subdir, 'header_nest'))
    
    # Get the domain definition
    species = getattr(controlvect.components, 'fluxes').parameters.attributes[0]
    tracer = getattr(getattr(controlvect.components, 'fluxes').parameters,
                     species)
    
    if tracer.hresol == 'regions':
        nbox = tracer.nregions
        
        # TODO: this will change once initial conditions are added (ciniopt)
        npvar = tracer.ndates * nbox
        ndvar = nbox
        nvar = npvar
        grad_o = np.zeros(nvar)
        
        ix1 = model.domain.ix1
        ix2 = model.domain.ix2
        iy1 = model.domain.iy1
        iy2 = model.domain.iy2
    else:
        raise Exception(
            "For FLEXPART, only hresol:regions are implemented in controlvect")
    
    # Loop through model periods and read model output
    self.missingperiod = False
    
    for di, df in zip(subsimu_dates[:-1], subsimu_dates[1:]):
        
        subdir = di.strftime("%Y%m")
        
        # Save to datastore for debugging purposes
        obs_ghg = np.empty(len(obsvect.datastore))
        obs_ghg[:] = np.nan
        obs_bkg = np.empty(len(obsvect.datastore))
        obs_bkg[:] = np.nan
        obs_sim = np.empty(len(obsvect.datastore))
        obs_sim[:] = np.nan
        obs_model = np.empty(len(obsvect.datastore))
        obs_model[:] = np.nan
        obs_check = np.empty(len(obsvect.datastore))
        obs_check[:] = np.nan
        obs_bkgerr = np.empty(len(obsvect.datastore))
        obs_bkgerr[:] = np.nan
        obs_err = np.empty(len(obsvect.datastore))
        obs_err[:] = np.nan
        
        # Loop over observation dates
        print("di, df", di, df, datetime.datetime.now())
        
        for obs_i, row in enumerate(obsvect.datastore.itertuples()):
            # For debugging
            obs_check[obs_i] = obs_i
            
            station = row.station
            
            runsubdir_nest = os.path.join(
                model.run_dir_nest, station.upper(), subdir)
            runsubdir_glob = os.path.join(
                model.run_dir_glob, station.upper(), subdir)
            
            file_date = row.Index.strftime('%Y%m%d%H%M%S')
            
            # Read nested grids
            file_name = 'grid_time_nest_' + file_date + '_001'
            
            if not os.path.isfile(os.path.join(runsubdir_nest, file_name)):
                continue
            
            grid_nest, gtime, ngrid = model.utils.read.read_flexpart_grid(
                runsubdir_nest, file_name, fp_header_nest,
                numscale=model.numscale)

            grid_nest *= model.coeff * model.mmair / model.molarmass
            
            # Read global footprints
            # TODO read correction factor dry air
            file_name = 'grid_time_' + file_date + '_001'
            
            if not os.path.isfile(os.path.join(runsubdir_glob, file_name)):
                continue
            
            grid_glob, gtime_glob, ngrid_glob = \
                model.utils.read.read_flexpart_grid(
                runsubdir_glob, file_name, fp_header_glob,
                numscale=model.numscale)
            
            grid_glob *= model.coeff * model.mmair / model.molarmass
            
            # Array for global transport background
            hbkg = np.sum(grid_glob[:, :, 0:ngrid - 1], axis=2)
            hbkg[ix1:ix2, iy1:iy2] = 0.0
            
            # Index to state vector
            ind = np.argmin(np.abs(tracer.dates[0::-1] - gtime_glob[0]))
            obs_bkg[obs_i] = np.sum(
                hbkg[:, :].T * controlvect.flx_bkg[ind, :, :])
            obs_bkgerr[obs_i] = obs_bkg[obs_i] * tracer.err
            
            # Transport for boxes in regions
            hbox = np.zeros((nbox * ngrid), np.float64)
            mask_reg = tracer.regions > 0
            for n in range(ngrid):
                inds = n * nbox + tracer.regions[mask_reg] - 1
                np.add.at(hbox, inds, grid_nest.T[n, mask_reg])
                
                # for jy in range(fp_header_nest.numy):
                #     for ix in range(fp_header_nest.numx):
                #         if tracer.regions[jy, ix] > 0:
                #             hbox[n*nbox + tracer.regions[jy, ix] - 1] += \
                #             grid_nest[ix, jy, n]
            
            # TODO: account for possible difference nested grid/inversion domain
            hnest = grid_nest
            
            istate = np.zeros(ngrid, dtype=int) - 1
            
            # Calculate indices to state vector
            for i, j in enumerate(gtime):
                if j > df:
                    istate[i] = -1
                else:
                    # Discard 1st tracer date in the comparison
                    mask = j - tracer.dates[1::] <= datetime.timedelta(0)
                    istate[i] = int(np.argmax(mask))
            
            if np.max(istate) < 0:
                continue
            
            # ntstep: number of inversion intervals covered by the footprints
            ntstep = int(
                np.max(istate[istate > -1]) - np.min(istate[istate > -1]) + 1)
            
            hx = np.zeros(ntstep * nbox)
            px = np.zeros(ntstep * nbox)
            
            for i in range(ngrid):
                if istate[i] == -1:
                    continue
                
                ns = istate[i] - np.min(istate[istate > -1]) + 1
                
                px[(ns - 1) * nbox:ns * nbox] = \
                    controlvect.x[istate[i] * nbox:(istate[i] + 1) * nbox]
                hx[(ns - 1) * nbox:ns * nbox] += hbox[i * nbox:(i + 1) * nbox]
            
            obs_model[obs_i] = np.dot(hx, px)
            
            # Change in mixing ratio from best guess estimate
            obs_ghg[obs_i] = 0.
            for i, itim in enumerate(gtime):
                ind = np.argmin(np.abs(tracer.dates[0::-1] - itim))
                obs_ghg[obs_i] += np.sum(
                    hnest.T[i, :, :] * controlvect.flxall[ind, :, :])
            
            if getattr(tracer, 'offsets', False):
                # Optimize offsets
                obs_sim[obs_i] = obs_model[obs_i] \
                                 + obs_ghg[obs_i] + obs_bkg[obs_i]
            else:
                # Optimize fluxes
                obs_sim[obs_i] = obs_model[obs_i] + obs_bkg[obs_i]
            
            # calculate gradient in observation space
            # Jo'(p) = H^TR^-1(H(p) - y)
            # calculate as: Jo'(p) = sum( H_i^T*ydel_i*R_i )
            
            # Contribution to observation gradient from obs_i
            departure = obs_sim[obs_i] - row.obs
            istate_uni = np.unique(istate).astype(int)
            
            for n in range(ntstep):
                if istate_uni[n] > -1:
                    grad_o[istate_uni[n] * ndvar:
                           (istate_uni[n] + 1) * ndvar] += \
                        hx[n * ndvar:(n + 1) * ndvar] \
                        * departure / (row.obserror ** 2
                                       + obs_bkgerr[obs_i] ** 2)
            
            print(obs_i, row.obs, obs_sim[obs_i])
            
        obsvect.dx = grad_o
        
        # Add the different components to datastore
        obsvect.datastore['obs_bkgerr'] = obs_bkgerr
        obsvect.datastore['sim'] = obs_sim
        obsvect.datastore['obs_ghg'] = obs_ghg
        obsvect.datastore['obs_bkg'] = obs_bkg
        obsvect.datastore['obs_model'] = obs_model
        obsvect.datastore['obs_sim'] = obs_sim
        obsvect.datastore['obs_check'] = obs_check
        obsvect.datastore['obs_err'] = np.sqrt(
            obsvect.datastore['obs_bkgerr'] ** 2 +
            obsvect.datastore['obserror'] ** 2)
        
        # Save grad_o for inspection 
        rundir = "{}/obsoperator".format(workdir)
        file_grado = '{}/grad_o_{}.txt'.format(rundir, run_id)
        np.savetxt(file_grado, grad_o, fmt='%.8e')
        
        model.output_read = True
        
        created = os.path.isdir(runsubdir_nest)
        
        # If the sub-directory was already created,
        # the observation operator considers that the sub-simulation
        # was already properly run, thus passing to next sub-periods
        # Compute the sub-simulation anyway if some previous periods
        # were missing (as stored in self.missingperiod)
        do_simu = (created
                   or not getattr(self, 'autorestart', False)
                   or self.missingperiod) and not read_only
        self.missingperiod = do_simu
        
        # Some verbose
        # verbose("Running {} for the period {} to {}"
        #         .format(model.plugin.name, di, df))
        # verbose("Running mode: {}".format(mode))
        # verbose("Sub-directory: {}".format(runsubdir_nest))
        # print "do_simu", do_simu
        # print "read_only", read_only
        
        # Prepare observations for the model
        if not read_only:
            model.dataobs = obsvect.obsvect2native(di, df, mode, runsubdir_nest,
                                                   workdir)
        
        # If only initializing inputs, continue to next sub-period
        if getattr(self, 'onlyinit', False):
            continue
        
        model.chain = min(di, df)
    
    # If only initializing inputs, exit
    if getattr(self, 'onlyinit', False):
        verbose("The run was correctly initialized")
        return
    
    # Re-initalizing the chain argument
    if hasattr(model, 'chain'):
        del model.chain
    
    if mode in ['fwd']:
        controlvect.dump(
            "{}/controlvect_{}".format(rundir, run_id),
            to_netcdf=getattr(controlvect, 'save_out_netcdf', True),
            dir_netcdf='{}/controlvect/'.format(rundir), run_id=run_id)
    
    rundir = "{}/obsoperator/{}".format(workdir, mode)
    path.init_dir(rundir)
    dump_type = obsvect.dump_type
    dump_datastore(obsvect.datastore,
                   file_monit='{}/monitor_{}_.{}'.format(rundir, run_id,
                                                         dump_type),
                   mode='w', dump_type=dump_type,
                   col2dump=['obs_ghg', 'obs_bkg', 'obs_model', 'obs_sim',
                             'obs_check',
                             'obs_bkgerr', 'obs_err', 'obs_hx'])
    
    # Returning obsvect to the simulator
    if mode is 'fwd':
        return obsvect
