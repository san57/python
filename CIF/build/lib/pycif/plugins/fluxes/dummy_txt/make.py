from __future__ import division

import numpy as np
import xarray as xr
from PIL import Image, ImageDraw, ImageFont
from builtins import map
from builtins import zip


def make(self, flx_file, flx_txt):

    # Split the text
    splits = split_text(flx_txt)

    # Domain parameters
    dom_xsize = self.domain.nlon
    dom_ysize = self.domain.nlat

    # If no text return zeros
    if flx_txt.strip() == "":
        return np.zeros((dom_ysize, dom_xsize))

    # Loop over all possible combinations of newlines
    file_font = getattr(
        self, "file_font", "/usr/share/fonts/dejavu/DejaVuSans-Bold.ttf"
    )

    areas = []
    text_widths = []
    text_heights = []
    valid_splits = []
    for split in splits:

        width = max(list(map(len, split)))
        size = min(dom_xsize // width, dom_ysize // (1 + len(split)))
        size += 3 * size // 4

        pil_font = ImageFont.truetype(file_font, size=size, encoding="unic")

        area = 0
        widths = []
        heights = []
        for s in split:
            w, h = pil_font.getsize(s)
            widths.append(w)
            heights.append(h)
            area += w * h

        valid_splits.append(split)
        areas.append(area)
        text_widths.append(widths)
        text_heights.append(heights)

    # create a blank canvas with extra space between lines
    canvas = Image.new("RGB", [dom_xsize, dom_ysize], (255, 255, 255))
    draw = ImageDraw.Draw(canvas)

    # Get the one with maximum occupation of space with minimum number of lines
    imax = np.where(np.array(areas) >= 0.9 * np.max(areas))[0]
    stds = [np.std(list(map(len, valid_splits[i]))) for i in imax]
    iimax = np.where(np.array(stds) == np.min(stds))[0]
    imax = [imax[i] for i in iimax]
    lens = [len(valid_splits[i]) for i in iimax]
    iimax = np.where(np.array(lens) == np.min(lens))[0]
    imax = [imax[i] for i in iimax]

    imax = imax[0]

    text_width = text_widths[imax]
    text_height = text_heights[imax]
    split = valid_splits[imax]

    # draw the text onto the canvas
    offset = (
        (dom_xsize - max(text_width)) // 2,
        (dom_ysize - sum(text_height)) // 2,
    )
    white = "#000000"

    width = offset[0] + max(text_width) // 2

    y_text = offset[1]
    for s, w, h in zip(split, text_width, text_height):
        draw.text((width - (w // 2), y_text), s, font=pil_font, fill=white)
        y_text += h

    # Convert the canvas into an array with values in [0, 1]
    flx = (255 - np.asarray(canvas)) / 255.0
    flx = xr.DataArray(
        flx.mean(axis=2)[np.newaxis, np.newaxis, ...],
        dims=("time", "lev", "lat", "lon"),
    )

    return flx


def split_text(txt):
    """Return all possible splits for a given text

    Args:
        txt (str): the text to split

    Returns:
        list(str)

    Notes: This function is not optimized and can take a long time for string
    of more than a few words

    """
    # If empty text
    if txt == "":
        return [txt]

    split = txt.split()

    # Recursively loop over substring of the text
    if len(split) == 1:
        return [[txt]]
    elif len(split) == 2:
        return [[txt], split]
    else:
        prev = [[split[0]] + b for b in split_text(" ".join(txt.split()[1:]))]
        post = [
            b + [split[-1]] for b in split_text(" ".join(txt.split()[:-1]))
        ]

        return [[txt]] + prev + post
