import numpy as np
import pandas as pd
import xarray as xr


def reindex(array, levels={}, method="ffill"):
    """

    :param array: input xarray to re-index
    :param levels: levels to re-index, with corresponding values
    :param method: method to use in the reindex
    :return:
    """

    var_in = array

    # Loop over levels to re-index
    for lev in levels:
        index_in = var_in[lev].data
        index_out = levels[lev]
        if type(index_out) == xr.core.dataarray.DataArray:
            index_out = index_out.data
        
        index_tmp = np.unique(np.append(index_in, index_out))
        df_index = pd.DataFrame(np.arange(len(index_in)), index=index_in)
        df_index = df_index.reindex(index_tmp, method="ffill")

        idim = var_in.dims.index(lev)
        shape_out = [
            s if k != idim else len(index_out)
            for k, s in enumerate(var_in.shape)
        ]

        var_out = np.empty(shape_out)
        slice_out = [
            slice(None) if k != idim else range(len(index_out))
            for k, s in enumerate(var_in.shape)
        ]
        slice_in = [
            slice(None) if k != idim else df_index.loc[index_out, 0].values
            for k, s in enumerate(var_in.shape)
        ]
        
        var_out[slice_out] = var_in.data[slice_in]
        
        coords_out = {
            coord_lev:
                index_out if coord_lev == lev
                else var_in[coord_lev].data
            for coord_lev in var_in.coords
        }

        var_in = xr.DataArray(var_out, coords=coords_out, dims=var_in.dims)

    return var_in
