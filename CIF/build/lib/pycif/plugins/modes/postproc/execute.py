import datetime
import os

import numpy as np

from pycif.utils.datastores.dump import read_datastore


def execute(self, **kwargs):
    """Performs a variational inversion given a minimizer method and a
    simulator (i.e. a function to minimize and its gradient)

    Args:
        self (Plugin): definition of the mode set-up

    """

    # Working directory
    workdir = self.workdir

    # Control vector
    statevect = self.statevect

    # Observation operator
    obsoper = self.obsoperator

    # Simulation window
    datei = self.datei
    datef = self.datef

    # Loads a monitor pickle
    if hasattr(self, "file_monitor"):
        obsvect = obsoper.obsvect
        obsvect.datastore = read_datastore(self.file_monitor)

    # Loads a control vector pickle
    if hasattr(self, "statevect_file"):
        statevect.load(self.statevect_file)

    # Fluxes to maps
    flx_maps = {}
    flx_maps_bg = {}
    for trcr in statevect.components.fluxes.parameters.attributes:
        tracer = getattr(statevect.components.fluxes.parameters, trcr)

        dates = tracer.dates[:-1]

        # Fetching the correct slice in the control vector
        tmp = statevect.x[tracer.xpointer: tracer.xpointer + tracer.dim]

        # Projecting to a map
        flx_map = statevect.utils.scale2map(
            tmp, tracer, dates, statevect.domain
        )

        # Changing time format for saving in netcdf
        flx_map["time"] = flx_map["time"].astype(datetime.datetime)

        flx_maps[trcr] = flx_map

        tmp = statevect.xb[tracer.xpointer: tracer.xpointer + tracer.dim]
        flx_map = statevect.utils.scale2map(
            tmp, tracer, dates, statevect.domain
        )
        flx_maps_bg[trcr] = flx_map

    # Save output as NetCDF
    flx_maps[trcr].to_netcdf(
        "{}/statevect_x.nc".format(os.path.dirname(self.statevect_file))
    )
    flx_maps_bg[trcr].to_netcdf(
        "{}/statevect_xb.nc".format(os.path.dirname(self.statevect_file))
    )

    ratio = flx_maps[trcr] - flx_maps_bg[trcr]
    ratio = ratio.where((ratio > -np.inf) * (ratio < np.inf))
    trcr = "CH4"
    ratio.to_netcdf(
        "{}/ratio_x-xb.nc".format(os.path.dirname(self.statevect_file))
    )

    # # Plot global average
    # nseconds = [d.total_seconds() for d in dates[1:] - dates[:-1]] + [5 *
    # 86400]
    # nseconds = 3600 * 24 * 365
    # x = (flx_maps['CH4'] * controlvect.domain.areas.T).sum(axis=(1,2))
    # xb = (flx_maps_bg['CH4'] * controlvect.domain.areas.T).sum(axis=(1,2))
    #
    # plt.plot(dates, x, c='r', label='posterior')
    # plt.plot(dates, xb, c='b', label='prior')
    # plt.legend()
    # plt.savefig('{}/figures/TotalBudget.pdf'.format(self.workdir),
    #             bbox='tight')
    # plt.close()
    #
    # # Save Control vector to NetCDF
    # dataset = xr.Dataset(data_vars={'CH4_bg': flx_maps_bg['CH4'],
    #                                 'CH4_post': flx_maps['CH4']})
    # dataset['time'] = list(dataset['time'].values)
    # dataset.to_netcdf('{}/outputs/controlvect_CH4_0023.nc'.format(
    # self.workdir))
