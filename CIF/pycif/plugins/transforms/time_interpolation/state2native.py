import numpy as np
import pandas as pd
import xarray as xr

try:
    import cPickle as pickle
except ImportError:
    import pickle


def state2native(
    transform,
    data,
    mapper,
    mod_input,
    di,
    df,
    mode,
    runsubdir,
    workdir,
    trans_mode,
    **kwargs
):

    ddi = min(di, df)

    xmod = data.datastore

    # Re-order input dates if not already done
    for trid in mapper["inputs"]:
        var_in = xmod[trid]["spec"]
        target_dates = mapper["outputs"][trid]["input_dates"][ddi]
        ntimes, nlev, nlat, nlon = var_in.shape
        time = np.sort(
            [pd.Timestamp(t).to_pydatetime() for t in var_in.time.data]
        )
        var_out = var_in.reindex({"time": time})

        target_dates_tmp = list(set(list(target_dates) + list(time)))
        target_dates_tmp.sort()

        if transform.method == "linear":
            df_index = pd.DataFrame(range(len(time)), index=time)
            df_index = (
                df_index.reindex(target_dates_tmp)
                .interpolate(method="index")
                .reindex(target_dates)
            )

            var_out = np.empty((len(target_dates), nlev, nlat, nlon))
            for k, dd in enumerate(target_dates):
                ind = df_index.iloc[k, 0]
                dmin = np.floor(ind).astype(int)
                wgt = ind - dmin
                var_out[k] = (
                    var_in[dmin] * (1 - wgt)
                    + var_in[min(dmin + 1, len(time) - 1)] * wgt
                )

        xmod[trid]["spec"] = xr.DataArray(
            var_out,
            coords={"time": target_dates},
            dims=("time", "lev", "lat", "lon"),
        )

    data.datastore = xmod

    return data
