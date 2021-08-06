#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""check sub-module

This module handles logging options and verbose levels.

"""
from __future__ import print_function

import logging
import os

from pycif.utils import path
from .coloredlog import ColorFormatter


def init_log(logfile, workdir, loglevel):
    """Initializes the log file for verbose outputs.

    Args:
        logfile (str): log file name
        workdir (str): directory where to create the logfile
        loglevel (int): level of verbosity.
                        2 for debug level, 1 for standard outputs

    Returns:
        (str, str) : (full path to the log file,
                      absolute path to the working directory)

    Notes: Beware that the function overwrites any existing log file.
    """

    # Turning the path to absolute and creating the directory
    workdir, _ = path.init_dir(workdir)
    if not os.path.isabs(logfile):
        logfile = "{}/{}".format(workdir, logfile)

    # Beware that the log_file is over-writen anyway
    # Flushing if exists
    open(logfile, "w").close()

    # Transform pycif verbose level to logging verbose level
    level = logging.DEBUG
    if loglevel == 2:
        level = logging.INFO
    elif loglevel > 2:
        level = logging.WARNING

    # Set up colored log
    stream_handler = logging.StreamHandler()
    formatter = ColorFormatter(fmt="#(level)%(message)s")
    stream_handler.setFormatter(formatter)

    file_handler = logging.FileHandler(logfile)
    formatter = ColorFormatter(
        fmt="%(asctime)s: %(message)s", datefmt="%Y-%m-%d %H:%m:%S"
    )
    file_handler.setFormatter(formatter)

    logger = logging.getLogger("")
    logger.setLevel(level)

    logger.addHandler(stream_handler)
    logger.addHandler(file_handler)

    return logfile, workdir


def debug(entry, logfile=None):
    return verbose(entry, logfile, verbose_level=1)


def info(entry, logfile=None):
    return verbose(entry, logfile, verbose_level=2)


def warning(entry, logfile=None):
    return verbose(entry, logfile, verbose_level=3)


def verbose(entry, logfile=None, verbose_level=2):
    """Prints out a log entry to the log_file

    Args:
        entry (string): entry to print
        logfile (file path): path to the log file
        verbose_level (int): level of verbosity

    Returns:
        None

    """
    verbose_threshold = logging.getLogger().level / 10

    if verbose_level < verbose_threshold:
        return

    if not isinstance(entry, list):
        entry = [entry]

    for ln in entry:
        if logfile is not None:
            with open(logfile, "a") as f:
                f.write(ln + "\n")

        else:
            if verbose_level == 1:
                logging.debug(ln)
            elif verbose_level == 2:
                logging.info(ln)
            if verbose_level > 3:
                logging.warning(ln)


#
# def check_memory(logfile=None):
#     info("Current memory usage: {} Mb"
#             .format(psutil.Process(os.getpid()).memory_info()[0]
#                     / float(2 ** 20)),
#             logfile)
