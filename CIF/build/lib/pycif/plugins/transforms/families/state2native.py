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

    keys = ["spec"]
    if mode == "tl":
        keys.append("incr")

    trid_out = mapper["outputs"].keys()[0]

    xmod[trid_out] = {k: 0 for k in keys}

    for trid in mapper["inputs"]:
        for k in keys:
            xmod[trid_out][k] += xmod[trid][k]

        for k in xmod[trid]:
            if k not in keys:
                xmod[trid_out][k] = xmod[trid][k]

        del xmod[trid]

    data.datastore = xmod

    return data
