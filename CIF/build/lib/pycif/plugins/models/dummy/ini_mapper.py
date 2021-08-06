def ini_mapper(model, transform_type, inputs={}, outputs={}, backup_comps={}):

    default_dict = {
        "input_dates": model.input_dates,
        "force_read": True,
        "force_dump": True,
    }
    dict_surface = dict(default_dict, **{"domain": model.domain})

    # Executable
    mapper = {
        "inputs": {("meteo", ""): dict_surface, ("param", ""): dict_surface},
        "outputs": {},
    }

    emis = {
        ("fluxes", s): dict_surface
        for s in model.chemistry.acspecies.attributes
    }
    mapper["inputs"].update(emis)

    return mapper
