def rinvprod(obsvect, dy):

    # At the moment, this is valid for diagonal observation matrices only
    # TODO: Include these lines into rinvprod.py, with option for
    # non-diagonal matrices, eventually
    return dy / obsvect.datastore["obserror"] ** 2
