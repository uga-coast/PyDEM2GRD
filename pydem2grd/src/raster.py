"""Raster coordinate and data helpers backed by Rasterio."""

from pathlib import Path

import rasterio
from rasterio.transform import rowcol, xy


def open_raster(path):
    """Open a raster path with Rasterio."""
    return rasterio.open(Path(path))


def get_rastvalue(col, row, rdata):
    """Return the first-band value at a zero-based column and row."""
    return rdata.read(1)[row, col]


def get_boundingbox(rdata):
    """Return raster bounds as ``[xmin, ymin, xmax, ymax]``."""
    bounds = rdata.bounds
    return [bounds.left, bounds.bottom, bounds.right, bounds.top]


def get_numrowcol(rdata):
    """Return the raster's column and row counts."""
    return rdata.width, rdata.height


def isinraster(x, y, rdata):
    """Return whether a coordinate lies strictly inside the raster bounds."""
    xmin, ymin, xmax, ymax = get_boundingbox(rdata)
    return xmin < x < xmax and ymin < y < ymax


def get_rastersize(rdata):
    """Return the absolute horizontal pixel resolution."""
    return abs(float(rdata.res[0]))


def coord2pixel(x, y, rdata):
    """Convert map coordinates to a zero-based ``(column, row)`` pair."""
    if not isinraster(x, y, rdata):
        return -1, -1
    row, col = rowcol(rdata.transform, x, y)
    return int(col), int(row)


def pixel2coord(col, row, rdata):
    """Return map coordinates for the upper-left corner of a raster cell."""
    xcoord, ycoord = xy(rdata.transform, row, col, offset="ul")
    return float(xcoord), float(ycoord)


def read_band_as_array(rdata, masked=False):
    """Read the first raster band."""
    return rdata.read(1, masked=masked)
