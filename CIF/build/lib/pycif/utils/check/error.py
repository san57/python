#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""check sub-module

This module handles logging options and verbose levels.

"""
from __future__ import print_function

import sys
import traceback
from functools import wraps

from . import warning


def handle_except(behave="exit"):
    """Decorates functions to handle exceptions in a generic way
    """

    def decorate(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            try:
                return func(*args, **kwargs)
            except Exception as e:
                bhv = kwargs.pop("behave", behave)
                error_handler(e, func, bhv, locargs=args, lockwargs=kwargs)

        return wrapper

    return decorate


def error_handler(
    error,
    fct,
    behave,
    locargs=[],
    lockwargs={},
    color="\033[91m",
    end="\033[0m",
    bold="\033[1m",
    *args,
    **kwargs
):
    """Attached to handle_except to define error messages.
    """

    if error.__class__ == IOError:
        towrite = "Couldn't open the file{}".format(error.message)
    elif error.__class__ == KeyError:
        towrite = "KeyValue {} is not known".format(error.message)
    else:
        towrite = error.message

    warning(bold + color + towrite + end)
    print(traceback.format_exc())

    if behave == "exit":
        sys.exit()
    elif behave == "pass":
        pass
    elif behave == "continue":
        return
    elif behave == "debug":
        # pdb.post_mortem()
        import code

        code.interact(local=dict(globals(), **locals()))
