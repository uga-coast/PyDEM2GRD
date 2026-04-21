from __future__ import annotations

from pathlib import Path
import numpy as np
import rasterio
from rasterio.transform import rowcol, xy


def open_raster(path: str | Path):
    return rasterio.open(path)


def get_boundingbox(rdata) -> list[float]:
    bounds = rdata.bounds
    return [bounds.left, bounds.bottom, bounds.right, bounds.top]


def get_numrowcol(rdata) -> tuple[int, int]:
    return rdata.width, rdata.height


def isinraster(x: float, y: float, rdata) -> bool:
    bbox = get_boundingbox(rdata)
    return bbox[0] < x < bbox[2] and bbox[1] < y < bbox[3]


def get_rastersize(rdata) -> float:
    return float(abs(rdata.transform.a))


def coord2pixel(x: float, y: float, rdata) -> tuple[int, int]:
    if not isinraster(x, y, rdata):
        return -1, -1
    row, col = rowcol(rdata.transform, x, y)
    return int(col), int(row)


def pixel2coord(col: int, row: int, rdata) -> tuple[float, float]:
    x, y = xy(rdata.transform, row, col, offset="center")
    return float(x), float(y)


def read_band_as_array(rdata) -> np.ndarray:
    return rdata.read(1)
