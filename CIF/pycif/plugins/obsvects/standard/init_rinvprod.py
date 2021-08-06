def init_rinvprod(obsvect, measurements, **kwargs):

    # Update observational errors by including transport errors
    # At the moment, transport error is specified by the user and considered
    # uniform
    obsvect.datastore.loc[
        obsvect.datastore["obserror"] <= 0, "obserror"
    ] = obsvect.datastore["obserror"].mean()
