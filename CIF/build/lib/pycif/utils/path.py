#!/usr/bin/env python
# -*- coding: utf-8 -*-

import errno
import glob
import os
import shutil


def init_dir(path):
    """Equivalent to "mkdir -p"

    Args:
        path (str): path to the directory to create

    Returns:
        True if created directory correctly, False if directory already existed

    """

    path = os.path.abspath(path)

    try:
        os.makedirs(path)
        return path, True

    # If already exists, pass; else, raise exception
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            return path, False
        else:
            raise


def link(source, target):
    """Equivalent to "ln -sf". Force overwriting if the target exists.

    Args:
        source (str): path to the source
        target (str): path to the target

    """

    # If the same do nothin
    if source == target:
        return

    # Unlinks if existing
    try:
        os.unlink(target)

    except OSError as e:
        if e.errno == errno.ENOENT:
            pass
        else:
            raise e

    # Raise error if origin does not exist
    if not os.path.exists(source):
        raise IOError(
            "Trying to link from a non-existing file: {}".format(source)
        )

    # Otherwise do the link
    os.symlink(source, target)


def remove(pattern):
    """Removes files and directories following a given pattern.
    Equivalent to 'rm -rf' """

    toremove = glob.glob(pattern)

    for path in toremove:
        if os.path.isfile(path):
            os.remove(path)
        elif os.path.isdir(path):
            shutil.rmtree(path)


def copyfromlink(fileout):
    # If file already exists but is link, first copy info from it
    if os.path.islink(fileout):
        orig_link = os.readlink(fileout)
        os.remove(fileout)
        shutil.copy(orig_link, fileout)
