def ini_mapper(model, transform_type, inputs={}, outputs={}, backup_comps={}):

    default_dict = {"input_dates": model.input_dates, "force_dump": True}
    dict_surface = dict(default_dict, **{"domain": model.domain})
    dict_ini = dict(
        dict_surface, **{"input_dates": {model.datei: [model.datei]}}
    )

    # Executable
    mapper = {
        "inputs": {
            ("exe", ""): default_dict,
            ("def", ""): default_dict,
            ("traj", ""): default_dict,
            ("meteo", ""): default_dict,
            ("chem_fields", ""): default_dict,
        },
        "outputs": {},
    }

    emis = {
        ("fluxes", s): dict_surface
        for s in model.chemistry.emis_species.attributes
    }
    inicond = {
        ("inicond", s): dict_ini for s in model.chemistry.acspecies.attributes
    }
    prescrcond = {
        ("prescrconcs", s): dict_ini
        for s in model.chemistry.prescrconcs.attributes
    }

    if hasattr(model.chemistry, "prodloss3d"):
        prodloss3d = {
            ("prodloss3d", s): dict_ini
            for s in model.chemistry.prodloss3d.attributes
        }
    else:
        prodloss3d = {}

    mapper["inputs"].update(
        dict(emis, **dict(inicond, **dict(prescrcond, **prodloss3d)))
    )

    # Accepts backup components instead of reference ones
    if hasattr(model, "backup_comps"):
        backup_comps.update(model.backup_comps)

    return mapper
