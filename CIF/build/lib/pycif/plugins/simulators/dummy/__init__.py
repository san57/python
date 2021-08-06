from .simul import simul

requirements = {
    "statevect": {"any": True, "empty": False},
    "obsvect": {"any": True, "empty": False},
    "obsoperator": {
        "any": True,
        "empty": True,
        "name": "standard",
        "version": "std",
    },
}
