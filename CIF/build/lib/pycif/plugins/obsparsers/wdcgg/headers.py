# -*- coding: utf-8 -*-

from __future__ import absolute_import
from __future__ import print_function

from pycif.utils.check import info
from .utils import remap_extract, find_header


def get_header(obs_file, maxlen):
    """Extract the header from a WDCGG File

    Args:
        obs_file (str): path to input file
        maxlen (int): abort after this amount of lines when reading header.
                      Default 300

    Returns:
        List[str]: List with all Lines of the Header

    """
    with open(obs_file, "r") as input_file:
        lines = []
        nheader = 0

        # Accepts formats with 'CXX' terminating the header
        # or with the number of HEADER LINES specified explicitly in the header
        for line in input_file:
            lines.append(line.strip())
            if "CXX" in line:
                return lines

            if "HEADER LINES" in line:
                nheader = int(line.split(":")[1])

            if len(lines) > maxlen:
                break

        if not lines:
            return []

        # if the number of line was not found, tries to find it
        # on the first line of the document
        if nheader == 0:
            try:
                nheader = int(lines[0].split()[1])

            except BaseException:
                raise ValueError(
                    "More than {} Header Lines in WDCGG File. "
                    "Is it a WDCGG file?".format(maxlen)
                )

        return lines[:nheader]


def parse_header(
    header, spec, list_extract, default_unit="ppm", default_tz="utc"
):
    """Extract information from the header

    Args:
        header (list[str]): extracted header
        spec (str): species to extract
        list_extract (list[str]): list of parameters to return
                                  'flag' to extract flag
                                  'error' to extract observation error
                                  any other parameter appearing in the columns
        default_unit (str): default unit generally used to report this species
        default_tz (str): default time zone for this file

    Returns:
        a 4-element tuple containing
            - names (list[str]): list of columns names to extract
            - columns (list[int]): list of column index to extract
            - date_ids (list[int]): list of column ids for date information
            - extra (dict): extra information contained in the header and not
                in the body of the file, e.g., altitude, coordinates, unit, etc.

    """
    # Minimize all characters to facilitate comparisons
    head = [s.lower() for s in header[-1].split()[1:]]

    # Parsing time information
    try:
        date_ids = [head.index("date")]

    except BaseException:
        info(header)
        info(head)
        raise ValueError(
            "Cant find a date in this WDCGG file. " "Please check format"
        )

    if "time" in head:
        date_ids.append(head.index("time"))

    # Getting other parameters using the utils.remap_extract function
    columns = []
    names = []

    extra = {}

    for id_extract in list_extract:
        try:
            # First look into columns names
            columns.append(head.index(remap_extract(id_extract)))
            names.append(id_extract.lower())

        except BaseException:
            try:
                # Some files have a name with CH4_Air instead of CH4
                columns.append(head.index(remap_extract(id_extract) + "_air"))
                names.append(id_extract.lower())

            except BaseException:
                try:
                    # Look into the header
                    id_value = find_header(id_extract, header)
                    extra[id_extract.lower()] = id_value

                except Exception as e:
                    # If cannot find,
                    # assume default values for unit and timezone
                    info("Cant extract " + id_extract)

                    if id_extract == "unit":
                        extra[id_extract] = default_unit

                    elif id_extract == "tz":
                        extra[id_extract] = default_tz

                    else:
                        extra[id_extract] = None

    return names, columns, date_ids, extra
