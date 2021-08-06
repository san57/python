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
    do_simu=True,
    **kwargs
):

    transform.model.run(runsubdir, mode, workdir, do_simu=do_simu, **kwargs)
    return data
