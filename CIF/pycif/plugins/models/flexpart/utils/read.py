""" Routines for reading FLEXPART output

"""

from __future__ import absolute_import
import os
import numpy as np
import glob

from scipy.io import FortranFile
import pandas as pd
import datetime

from . import flexpart_header
from . import mod_flexpart
# from pycif.utils.dates import j2d


def read_flexpart_dir(subdir, nested=True, **kwargs):
    """ Reads all FLEXPART grid_time files from a directory

    Returns:
        grid_fp: Array (nlon, nlat, ntime, nread) of footprints
        gtimes : Array (nread) of file dates
    
    """

    header_filename='header_nest'
    if not nested:
        header_filename='header'
        
    fp_header = flexpart_header.Flexpartheader()
    fp_header.read_header(os.path.join(subdir, header_filename))
    file_dates = os.path.join(subdir, 'dates')

    if nested:
        files_list = glob.glob('{}/grid_time_nest_*'.format(subdir))
    else:
        files_list = glob.glob('{}/grid_time_[0-9]*'.format(subdir))
        
    nread = len(files_list)

    grid = np.ndarray((fp_header.numx, fp_header.numy, fp_header.maxngrid))
    grids = np.ndarray((fp_header.numx, fp_header.numy, fp_header.maxngrid, nread))
    gtimes = np.ndarray((fp_header.maxngrid, nread))
    
    for ii, ifile in enumerate(files_list):
        grid_fp, ngrid, gtime = mod_flexpart.read_grid(
        ifile, file_dates, fp_header.numx, fp_header.numy, \
        fp_header.maxngrid, fp_header.xshift,  fp_header.ndgrid, fp_header.trajdays)

        grids[:,:,:,ii] = grid_fp
        gtimes[:, ii] = gtime

        
    return grids, gtimes


def read_flexpart_grid(subdir, file_name, fp_header, **kwargs):
    """ Reads FLEXPART grid_time file.
    Convert s.m3/kg to s.m2/kg
    

    Inputs:
        file_name - full path name to file name
        fp_header - a Flexpartheader object

    Returns:
        grid_fp: Array (nlon, nlat, ntime) of footprints
        gtime  : Array (maxngrid) of file dates
    
    """

    # Numeric scaling factor
    numscale = np.float(kwargs.get("numscale", 1.E12))
        
    file_dates = os.path.join(subdir, 'dates')
    path_file = os.path.join(subdir, file_name)

    scaleconc = 1.e12
    
    days = []
    times = []
    counts_i = []
    counts_r = []
    sparse_i = []
    sparse_r = []
    with FortranFile(path_file, 'r') as f:
        # Looping until the end of the binary file
        while True:
            try:
                yyyymmdd = f.read_ints('i4')[0].astype(str)
                hhmmss = '{0:06d}'.format(f.read_ints('i4')[0])
                i = f.read_ints('i4')
                spi = f.read_ints('i4')
                r = f.read_ints('i4')
                spr = f.read_reals('f4')
                
                days.append(yyyymmdd)
                times.append(hhmmss)
                counts_i.extend(i)
                sparse_i.extend(spi)
                counts_r.extend(r)
                sparse_r.extend(spr)
    
            except TypeError:
                break

    gtime_dt = [datetime.datetime.strptime('{}{}'.format(d, h), '%Y%m%d%H%M%S')
                for d, h in zip(days, times)][::-1]
    ngrid = len(gtime_dt)
    grid_fp = np.zeros((fp_header.numx, fp_header.numy, ngrid + 2))
    
    if len(sparse_r) > 0:
        sign = np.sign(sparse_r)
        cum_counts = np.unique(np.cumsum(counts_r) % sign.size)
        signchange = ((np.roll(sign, 1) - sign) != 0)
        signchange[cum_counts] = True
        
        inds_out = np.zeros((len(sparse_r)), dtype=np.int)
        inds_out[signchange] = sparse_i
        
        mask = inds_out == 0
        idx = np.where(~mask, np.arange(mask.size), 0)
        np.maximum.accumulate(idx, out=idx)
        inds_out = inds_out[idx] + np.arange(mask.size) \
                   - idx - fp_header.numx * fp_header.numy
        jy, jx = np.unravel_index(inds_out, (fp_header.numy, fp_header.numx))
        
        jt = np.zeros((len(sparse_r)), dtype=np.int)
        jt[cum_counts] = np.arange(ngrid)[np.array(counts_r) > 0]
        np.maximum.accumulate(jt, out=jt)
        jt = ngrid - jt - 1
        grid_fp[jx, jy, jt] = np.abs(sparse_r) * scaleconc
    
    # grid_fp, ngrid, gtime = mod_flexpart.read_grid(
    #     path_file, file_dates, fp_header.numx, fp_header.numy,
    #     fp_header.maxngrid, fp_header.xshift,  fp_header.ndgrid, fp_header.trajdays)

    # Convert from ppt to ppmv or ppbv

    # Convert s.m3/kg to s.m2/kg and apply numerical scaling
    grid_fp = grid_fp / (fp_header.outheight[0]*numscale)

#     # Convert grid times to datetime format
#     gtime_dt = []
#     for i in range(len(gtime)):
# #        gtime_dt[i] = flexpart_header.Flexpartheader.j2d(gtime[i])
#         if gtime[i] == 0.:
#             break
#
# #        gtime_dt.append(fp_header.j2d(gtime[i]))
#         gtime_dt.append(j2d(gtime[i]))
    
    return grid_fp, gtime_dt, ngrid


def get_spec(subdir, **kwargs):
    """ Get species name from simulation output

    """

    with open(os.path.join(subdir, 'header_txt'), 'r') as f:
        txt = f.readlines()

    spec = txt[20].split()[1]

    print("species", spec)
    return spec
