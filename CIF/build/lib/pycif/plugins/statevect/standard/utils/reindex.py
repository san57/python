import pandas as pd


def reindex(array, var, levels={}, method="ffill"):
    """

    :param array: input xarray to re-index
    :param var: field in array to reindex
    :param levels: levels to re-index, with corresponding values
    :param method: method to use in the reindex
    :return:
    """

    dataframe = array[var].to_dataframe(name=var)
    df_levs = list(dataframe.index.names)

    # Loop over levels to re-index
    for lev in levels:
        loc_levs = list(set(df_levs) - {lev})

        lev_data = (
            dataframe.unstack(level=loc_levs)
            .reindex(pd.Index(levels[lev], name=lev), method=method)
            .stack(level=loc_levs)
        )
        dataframe = lev_data

    # Get back to xarray type, re-order dimensions to fit input data
    dataframe = dataframe.to_xarray()[var].transpose(*df_levs)

    return dataframe
