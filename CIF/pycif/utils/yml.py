import os
import re
from collections import OrderedDict

import yaml


def ordered_load(
    stream, Loader=yaml.SafeLoader, object_pairs_hook=OrderedDict
):
    """A function from Stack-OverFlow to keep the order of lines in the Yaml.
    See: https://stackoverflow.com/questions/5121931/in-python-how-can-you-load-yaml-mappings-as-ordereddicts"""

    # Ordered loader
    class OrderedLoader(Loader):
        pass

    def construct_mapping(loader, node):
        loader.flatten_mapping(node)
        return object_pairs_hook(loader.construct_pairs(node))

    OrderedLoader.add_constructor(
        yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG, construct_mapping
    )

    # Resolving bash environment variables
    path_matcher = re.compile(r"(\$\{([^}^{]+)\}|~)")

    def path_constructor(loader, node):
        """Extract the matched value,
        expand env variable, and replace the match.
        If home tilde is used expanduser"""
        value = node.value
        match = path_matcher.match(value)
        group = match.group()
        if group == "~":
            return os.path.expanduser(value)

        env_var = group[2:-1]
        return os.environ.get(env_var) + value[match.end():]

    OrderedLoader.add_implicit_resolver("!path", path_matcher, None)
    OrderedLoader.add_constructor("!path", path_constructor)

    return yaml.load(stream, OrderedLoader)
