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

    xmod = data.datastore
    for trid in mapper["inputs"]:
        xmod[trid]["spec"] *= transform.scale

    data.datastore = xmod
    return data
