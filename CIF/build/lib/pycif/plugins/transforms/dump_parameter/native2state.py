from pycif.utils.check import info


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

    datastore = data.datastore
    for trid in mapper["inputs"]:
        input_type = trid[0]

        # If trid in datastore, dumps this one
        if trid in datastore:
            todump = [trid]

        # If input parameter is '',
        # dumps all available parameters of this component
        elif trid[1] == "":
            todump = [t for t in datastore if t[0] == trid[0]]

        # Otherwise check whether there is a component
        # encompassing all parameters (i.e., with '' as parameter)
        else:
            todump = [(trid[0], "")]

        # Create new data to extract
        data2dump = {t: datastore[t] for t in todump}

        data2dump = transform.model.outputs2native(
            data2dump, input_type, di, df, runsubdir, mode
        )

        for tr in data2dump:
            if tr in datastore:
                datastore[tr].update(data2dump[tr])

            else:
                info(
                    "{} was simulated by the model, "
                    "but could not be transferred to the control vector".format(
                        tr
                    )
                )

    return data
