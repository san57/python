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

    datastore = data.datastore
    for trid in mapper["inputs"]:
        input_type = trid[0]

        # If trid in datastore, dumps this one
        if trid in datastore:
            todump = [trid]

        # If input parameter is '',
        # dumps all available parameters of this component
        elif trid[1] == "":
            todump = [t for t in datastore if t[0] == input_type]

        # Otherwise check whether there is a component
        # encompassing all parameters (i.e., with '' as parameter)
        else:
            todump = [(input_type, "")]

        # Create new data to dump
        data2dump = transform.from_dict({})
        data2dump.datastore = {
            t: datastore[t] for t in todump if t in datastore
        }

        # If the model does not need to compute a simulation,
        # just skip this step
        if do_simu:
            transform.model.native2inputs(
                data2dump, input_type, di, df, runsubdir, mode
            )

        # Dumping data from the datastore
        for t in todump:
            if t not in datastore:
                continue

            if "spec" in datastore[t]:
                del datastore[t]["spec"]
            if "incr" in datastore[t]:
                del datastore[t]["incr"]

    return data
