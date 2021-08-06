import copy


def native2state(
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

    xmod = data.datastore
    trid_out = mapper["outputs"].keys()[0]

    for trid in mapper["inputs"]:
        xmod[trid] = {k: xmod[trid_out][k] for k in xmod[trid_out]}
        xmod[trid]["adj_out"] = copy.deepcopy(xmod[trid_out]["adj_out"])

    del xmod[trid_out]

    data.datastore = xmod

    return data
