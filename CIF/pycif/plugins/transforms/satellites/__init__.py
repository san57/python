from __future__ import absolute_import

from .native2obsvect import native2obsvect
from .obsvect2native import obsvect2native

requirements = {"model": {"any": True, "empty": False}}


def ini_mapper(
    transf, mapper,
):
    """Initialize mapper for satellites.
    Does not need to change the mapper as cropping time keep similar
    names and dimensions"""

    return mapper
