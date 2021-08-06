def ini_mapper(model, transform_type, inputs={}, outputs={}, backup_comps={}):

    default_dict = {"input_dates": model.input_dates, "force_dump": True}
    dict_surface = dict(default_dict, **{"domain": model.domain})
    dict_bound = dict(dict_surface, **{"is_lbc": True})
    dict_top = dict(dict_surface, **{"is_top": True})
    dict_ini = dict(
        dict_surface, **{"input_dates": {model.datei: [model.datei]}}
    )

    # Executable
    mapper = {
        "inputs": {
            ("exe", ""): default_dict,
            ("nml", ""): default_dict,
            ("meteo", ""): default_dict,
        },
        "outputs": {
            ("concs", s): {} for s in model.chemistry.acspecies.attributes
        },
    }

    emis = {
        ("fluxes", s): dict_surface
        for s in model.chemistry.emis_species.attributes
    }
    bioemis = {
        ("biofluxes", s): dict_surface
        for s in model.chemistry.bio_species.attributes
    }
    inicond = {
        ("inicond", s): dict_ini for s in model.chemistry.acspecies.attributes
    }
    lbc = {
        ("latcond", s): dict_bound
        for s in model.chemistry.acspecies.attributes
    }
    top = {
        ("topcond", s): dict_top for s in model.chemistry.acspecies.attributes
    }

    mapper["inputs"].update({**emis, **bioemis, **inicond, **lbc, **top})

    # Accepts backup components instead of reference ones
    backup_comps.update(model.backup_comps)

    return mapper
