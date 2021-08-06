import numpy as np

from pycif.utils.classes.setup import Setup


def add_default(
    transforms,
    yml_dict,
    position="last",
    index=0,
    init=False,
    mapper={},
    transform_type="state",
    inputs={},
    outputs={},
    backup_comps={},
):
    # Create a random name to identify the current transformation
    new_id = str(np.random.randint(0, 1e10, 1)[0])
    while hasattr(transforms, new_id):
        new_id = str(np.random.randint(0, 1e10, 1)[0])

    new_transf = Setup.from_dict({new_id: yml_dict})
    Setup.load_setup(new_transf)
    new_transf = getattr(new_transf, new_id)
    setattr(transforms, new_id, new_transf)

    if position == "last":
        transforms.attributes.append(new_id)

    elif position == "start":
        transforms.attributes.insert(0, new_id)

    elif position == "index":
        transforms.attributes.insert(index, new_id)

    # Initializes mapper
    if init:
        transf_mapper_loc = new_transf.ini_mapper(
            transform_type,
            inputs=inputs,
            outputs=outputs,
            backup_comps=backup_comps,
        )

        new_transf.mapper = transf_mapper_loc
        mapper[new_id] = transf_mapper_loc

        inputs.update(transf_mapper_loc["inputs"])
        outputs.update(transf_mapper_loc["outputs"])

    return new_transf, new_id
