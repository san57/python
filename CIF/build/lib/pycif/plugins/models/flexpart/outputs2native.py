import numpy as np
import datetime
import os
import calendar
import glob
import pandas as pd
import xarray as xr
from .utils import read

from pycif.utils.netcdf import save_nc
from pycif.utils.check import info, debug
from pycif.utils.datastores.empty import init_empty

import pickle


def outputs2native(self, data2dump, input_type,
                   di, df, runsubdir, mode='fwd', **kwargs):
    """Reads outputs to pycif objects.
       
       Does nothing for now as we instead read FLEXPART output 
       inside loop over observations in obsoper.py 

    """
    
    ddi = min(di, df)
    
    dataobs = self.dataobs
    nobs = len(dataobs)
    subdir = ddi.strftime("%Y%m")
    
    # Initialize header
    fp_header_glob = self.utils.flexpart_header.Flexpartheader()
    fp_header_glob.read_header(
        os.path.join(self.run_dir_glob,
                     dataobs.head(1)['station'][0].upper(),
                     subdir, 'header'))
    
    fp_header_nest = None
    if self.nested:
        fp_header_nest = self.utils.flexpart_header.Flexpartheader()
        fp_header_nest.read_header(
            os.path.join(self.run_dir_nest,
                         dataobs.head(1)['station'][0].upper(),
                         subdir, 'header_nest'))
    
    # Nest domain definition
    ix1 = self.domain.ix1
    ix2 = self.domain.ix2
    iy1 = self.domain.iy1
    iy2 = self.domain.iy2
    
    # Save to datastore for debugging purposes
    obs_ghg = np.nan * np.empty(nobs)
    obs_bkg = np.nan * np.empty(nobs)
    obs_sim = np.nan * np.empty(nobs)
    obs_model = np.nan * np.empty(nobs)
    obs_check = np.nan * np.empty(nobs)
    obs_bkgerr = np.nan * np.empty(nobs)
    obs_err = np.nan * np.empty(nobs)
    
    # Loop over observation dates
    info("di, df: {}, {}, {}".format(di, df, datetime.datetime.now()))
    
    # Initialize output sensitivity
    if mode == "adj":
        nlat, nlon = self.domain.zlat.shape
        datasensit = np.zeros((len(self.input_dates[ddi]), 1, nlat, nlon))
    
    for obs_i, row in enumerate(dataobs.itertuples()):
        station = row.station
        
        # Read nested grids
        runsubdir_nest = os.path.join(
            self.run_dir_nest, station.upper(), subdir)
        file_date = row.Index.strftime('%Y%m%d%H%M%S')
        file_name = 'grid_time_nest_{}_001'.format(file_date)
        
        debug("Reading {}".format(file_name))
        
        if not os.path.isfile(os.path.join(runsubdir_nest, file_name)):
            continue
        
        grid_nest, gtime, ngrid = self.utils.read.read_flexpart_grid(
            runsubdir_nest, file_name, fp_header_nest,
            numscale=self.numscale)
        
        # Read global footprints
        # TODO read correction factor dry air)
        runsubdir_glob = os.path.join(
            self.run_dir_glob, station.upper(), subdir)
        file_name = 'grid_time_{}_001'.format(file_date)

        if not os.path.isfile(os.path.join(runsubdir_glob, file_name)):
            continue

        grid_glob, gtime_glob, ngrid_glob = \
            self.utils.read.read_flexpart_grid(
                runsubdir_glob, file_name, fp_header_glob,
                numscale=self.numscale)
        
        # Conversion of footprints
        grid_nest *= self.coeff * self.mmair / self.molarmass
        grid_glob *= self.coeff * self.mmair / self.molarmass
        
        # Find time stamps in dataflx to compare with grid
        # TODO: deal with several species
        fluxes = {"sim": list(self.dataflx.values())[0]["spec"]}
        if mode == "tl":
            fluxes["sim_tl"] = list(self.dataflx.values())[0]["incr"]
        
        # Multiple footprints with fluxes for fwd and tl
        if mode in ["fwd", "tl"]:
            for data_id in fluxes:
                dataflx = fluxes[data_id]
                flx_dates = pd.DatetimeIndex(
                    dataflx.time.values).to_pydatetime()
                
                inds_flx = \
                    np.argmin(np.abs(np.array(gtime)[:, np.newaxis]
                                     - flx_dates[np.newaxis, :]), axis=1)
                nest_sim = (
                        grid_nest.T[:ngrid].reshape(ngrid, -1)
                        * dataflx[inds_flx, 0, :self.domain.zlon_in.size, 0])\
                    .sum()
    
                inds_flx_glob = \
                    np.argmin(np.abs(np.array(gtime_glob)[:, np.newaxis]
                                     - flx_dates[np.newaxis, :]), axis=1)
                glob_sim = (
                        np.delete(grid_glob.T[:ngrid_glob].reshape(ngrid_glob, -1),
                                  self.domain.raveled_indexes_glob, 1)
                        * dataflx[inds_flx_glob, 0,
                                  self.domain.zlon_in.size:, 0]).sum()
                
                # Filling simulation
                sim_col = dataobs.columns.get_indexer([data_id])
                dataobs.iloc[obs_i, sim_col] = nest_sim + glob_sim
        
        # Adjoint
        else:
            dataflx = list(self.dataflx.values())[0]["spec"]
            flx_dates = pd.DatetimeIndex(
                dataflx.time.values).to_pydatetime()
            
            # Nest sensitivity
            nest_sensit = grid_nest.T[:ngrid].reshape(ngrid, -1)
            inds_flx = \
                np.argmin(np.abs(np.array(gtime)[:, np.newaxis]
                                 - flx_dates[np.newaxis, :]), axis=1)
            
            zeros = np.zeros((inds_flx.size, self.domain.zlon_in.size),
                             dtype=np.int)
            np.add.at(
                datasensit,
                (inds_flx.reshape(-1, 1),
                 zeros,
                 np.arange(self.domain.zlon_in.size)[np.newaxis, :],
                 zeros),
                nest_sensit * row.obs_incr,
            )
            
            # Global sensitivity
            glob_sensit = \
                np.delete(grid_glob.T[:ngrid_glob].reshape(ngrid_glob, -1),
                          self.domain.raveled_indexes_glob, 1)
            inds_flx_glob = \
                np.argmin(np.abs(np.array(gtime_glob)[:, np.newaxis]
                                 - flx_dates[np.newaxis, :]), axis=1)

            zeros = np.zeros((inds_flx.size,
                              nlat - self.domain.zlon_in.size),
                             dtype=np.int)
            np.add.at(
                datasensit,
                (inds_flx_glob.reshape(-1, 1),
                 zeros,
                 np.arange(self.domain.zlon_in.size,
                           nlat)[np.newaxis, :],
                 zeros),
                glob_sensit * row.obs_incr,
            )
            
    if mode in ["tl", "fwd"]:
        return self.dataobs
    
    else:
        outkey = list(data2dump.keys())[0]
        data2dump[outkey]["adj_out"] = xr.DataArray(
                        datasensit,
                        coords={"time": self.input_dates[ddi]},
                        dims=("time", "lev", "lat", "lon"))
        return data2dump
        
            
            # # Read global footprints
            # # TODO read correction factor dry air
            # file_name = 'grid_time_' + file_date + '_001'
            #
            # if not os.path.isfile(os.path.join(runsubdir_glob, file_name)):
            #     continue
            #
            # grid_glob, gtime_glob, ngrid_glob = \
            #     self.utils.read.read_flexpart_grid(
            #         runsubdir_glob, file_name, fp_header_glob,
            #         numscale=self.numscale)
            #
            # grid_glob *= self.coeff * self.mmair / self.molarmass
            #
            # # Array for global transport background
            # hbkg = np.sum(grid_glob[:, :, 0:ngrid_glob - 1], axis=2)
            # hbkg[ix1:ix2, iy1:iy2] = 0.0
            #
            # print(__file__)
            # import code
            # code.interact(local=dict(locals(), **globals()))
            # # Index to state vector
            # ind = np.argmin(np.abs(tracer.dates[0::-1] - gtime_glob[0]))
            # obs_bkg[obs_i] = np.sum(
            #     hbkg[:, :].T * controlvect.flx_bkg[ind, :, :])
            # obs_bkgerr[obs_i] = obs_bkg[obs_i] * tracer.err
            #
            # # Transport for boxes in regions
            # hbox = np.zeros((nbox * ngrid), np.float64)
            # mask_reg = tracer.regions > 0
            # for n in range(ngrid):
            #     inds = n * nbox + tracer.regions[mask_reg] - 1
            #     np.add.at(hbox, inds, grid_nest.T[n, mask_reg])
            #
            #     # for jy in range(fp_header_nest.numy):
            #     #     for ix in range(fp_header_nest.numx):
            #     #         if tracer.regions[jy, ix] > 0:
            #     #             hbox[n*nbox + tracer.regions[jy, ix] - 1] += \
            #     #             grid_nest[ix, jy, n]
            #
            # # TODO: account for possible difference nested grid/inversion domain
            # hnest = grid_nest
            #
            # istate = np.zeros(ngrid, dtype=int) - 1
            #
            # # Calculate indices to state vector
            # for i, j in enumerate(gtime):
            #     if j > df:
            #         istate[i] = -1
            #     else:
            #         # Discard 1st tracer date in the comparison
            #         mask = j - tracer.dates[1::] <= datetime.timedelta(0)
            #         istate[i] = int(np.argmax(mask))
            #
            # if np.max(istate) < 0:
            #     continue
            #
            # # ntstep: number of inversion intervals covered by the footprints
            # ntstep = int(
            #     np.max(istate[istate > -1]) - np.min(istate[istate > -1]) + 1)
            #
            # hx = np.zeros(ntstep * nbox)
            # px = np.zeros(ntstep * nbox)
            #
            # for i in range(ngrid):
            #     if istate[i] == -1:
            #         continue
            #
            #     ns = istate[i] - np.min(istate[istate > -1]) + 1
            #
            #     px[(ns - 1) * nbox:ns * nbox] = \
            #         controlvect.x[istate[i] * nbox:(istate[i] + 1) * nbox]
            #     hx[(ns - 1) * nbox:ns * nbox] += hbox[i * nbox:(i + 1) * nbox]
            #
            # obs_model[obs_i] = np.dot(hx, px)
            #
            # # Change in mixing ratio from best guess estimate
            # obs_ghg[obs_i] = 0.
            # for i, itim in enumerate(gtime):
            #     ind = np.argmin(np.abs(tracer.dates[0::-1] - itim))
            #     obs_ghg[obs_i] += np.sum(
            #         hnest.T[i, :, :] * controlvect.flxall[ind, :, :])
            #
            # if getattr(tracer, 'offsets', False):
            #     # Optimize offsets
            #     obs_sim[obs_i] = obs_model[obs_i] \
            #                      + obs_ghg[obs_i] + obs_bkg[obs_i]
            # else:
            #     # Optimize fluxes
            #     obs_sim[obs_i] = obs_model[obs_i] + obs_bkg[obs_i]
            #
            # # calculate gradient in observation space
            # # Jo'(p) = H^TR^-1(H(p) - y)
            # # calculate as: Jo'(p) = sum( H_i^T*ydel_i*R_i )
            #
            # # Contribution to observation gradient from obs_i
            # departure = obs_sim[obs_i] - row.obs
            # istate_uni = np.unique(istate).astype(int)
            #
            # for n in range(ntstep):
            #     if istate_uni[n] > -1:
            #         grad_o[istate_uni[n] * ndvar:
            #                (istate_uni[n] + 1) * ndvar] += \
            #             hx[n * ndvar:(n + 1) * ndvar] \
            #             * departure / (row.obserror ** 2
            #                            + obs_bkgerr[obs_i] ** 2)
    
