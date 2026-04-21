from __future__ import annotations

import math
import operator
from functools import reduce
from pathlib import Path
import numpy as np
import rasterio
from shapely.geometry import Point, Polygon, box, mapping
from rasterio.mask import mask

from .raster import (
    coord2pixel,
    get_boundingbox,
    get_numrowcol,
    get_rastersize,
    open_raster,
    pixel2coord,
    read_band_as_array,
)


def griddata(mesh, meshconn, xc, yc, boundaryNodes, raster, mfac, values, numvaluesgathered):
    mesh.size = mesh.computeMeshSize()

    with open_raster(raster) as data:
        bbox = get_boundingbox(data)
        bboxPoly = box(bbox[0], bbox[1], bbox[2], bbox[3])
        rastersize = get_rastersize(data)

    for i in range(mesh.numNodes()):
        bufr = 1.25 * mesh.size[i] * rastersize
        if (not bboxPoly.buffer(bufr).contains(Point(mesh.node(i).x(), mesh.node(i).y()))) or (
            mesh.node(i).z() > -999.0
        ):
            continue

        numElem = mesh.numElementsAroundNode(i)
        pointList = []
        for j in range(numElem):
            element = mesh.elementTable(mesh.node(i), j)
            pointList.append((xc[element.id() - 1], yc[element.id() - 1]))

        if mesh.node(i).id() in boundaryNodes:
            if not bboxPoly.contains(Point(mesh.node(i).x(), mesh.node(i).y())):
                continue
            pointList.append((mesh.node(i).x(), mesh.node(i).y()))

            if numElem == 1:
                element = mesh.elementTable(mesh.node(i), 0)
                for k in range(3):
                    if meshconn[element.id() - 1][k] == mesh.node(i).id():
                        continue
                    neigh = mesh.nodeIndexById(meshconn[element.id() - 1][k])
                    xmid = 0.5 * (mesh.node(i).x() + mesh.node(neigh).x())
                    ymid = 0.5 * (mesh.node(i).y() + mesh.node(neigh).y())
                    pointList.append((xmid, ymid))

        center = tuple(
            map(
                operator.truediv,
                reduce(lambda x, y: list(map(operator.add, x, y)), pointList),
                [len(pointList)] * 2,
            )
        )
        pointList = sorted(
            pointList,
            key=lambda coord: (-135 - math.degrees(math.atan2(*tuple(map(operator.sub, coord, center))[::-1]))) % 360,
        )
        vor = Polygon(pointList)
        geoms = [mapping(vor)]

        with rasterio.open(raster) as src:
            rbbox = Polygon(
                [
                    (src.bounds.left, src.bounds.bottom),
                    (src.bounds.left, src.bounds.top),
                    (src.bounds.right, src.bounds.top),
                    (src.bounds.right, src.bounds.bottom),
                ]
            )
            if not rbbox.intersects(vor):
                continue
            ndv = src.nodata
            out_image, _ = mask(src, geoms, crop=True)

        subset = np.asarray(out_image)
        if ndv is not None:
            subset = subset[subset > ndv]
        subset = subset * mfac
        numvaluesgathered[i] = numvaluesgathered[i] + subset.size
        values[i] = values[i] + np.sum(subset)

    return values, numvaluesgathered


def meshconnectivity(mesh):
    meshconn = mesh.connectivity()
    xc = np.zeros(mesh.numElements())
    yc = np.zeros(mesh.numElements())
    for i in range(mesh.numElements()):
        x1 = mesh.node(mesh.nodeIndexById(meshconn[i][0])).x()
        x2 = mesh.node(mesh.nodeIndexById(meshconn[i][1])).x()
        x3 = mesh.node(mesh.nodeIndexById(meshconn[i][2])).x()

        y1 = mesh.node(mesh.nodeIndexById(meshconn[i][0])).y()
        y2 = mesh.node(mesh.nodeIndexById(meshconn[i][1])).y()
        y3 = mesh.node(mesh.nodeIndexById(meshconn[i][2])).y()

        xc[i] = (x1 + x2 + x3) / 3.0
        yc[i] = (y1 + y2 + y3) / 3.0

    boundaryNodes = [int(node.id()) for node in mesh.boundaryNodes()]
    return meshconn, xc, yc, boundaryNodes


def gathervalues(mesh, raster, N, CA, mfac, values, numvaluesgathered):
    with open_raster(raster) as data:
        numcols, numrows = get_numrowcol(data)
        vals = read_band_as_array(data).astype(float)
        ndv = data.nodata
        vals = vals * mfac
        rastersize = get_rastersize(data)
        bbox = get_boundingbox(data)
        bboxPoly = box(bbox[0], bbox[1], bbox[2], bbox[3])

        for i in range(mesh.numNodes()):
            if mesh.node(i).z() > -999.0:
                continue

            if (
                mesh.node(i).x() < bbox[0] - (N[i] + 1 * rastersize)
                or mesh.node(i).x() > bbox[2] + (N[i] + 1 * rastersize)
                or mesh.node(i).y() < bbox[1] - (N[i] + 1 * rastersize)
                or mesh.node(i).y() > bbox[3] + (N[i] + 1 * rastersize)
            ):
                continue

            if numvaluesgathered[i] == CA[i]:
                continue

            col, row = coord2pixel(mesh.node(i).x(), mesh.node(i).y(), data)
            if row < 0 or col < 0:
                continue

            left = max(int(col - N[i]), 0)
            bottom = max(int(row + N[i]), 0)
            right = min(int(col + N[i]), numcols)
            top = max(int(row - N[i]), 0)

            xmin, ymin = pixel2coord(left, bottom, data)
            xmax, ymax = pixel2coord(right, top, data)
            stencilPoly = box(xmin, ymin, xmax, ymax)
            if stencilPoly.touches(bboxPoly):
                continue

            subset = vals[top:bottom, left:right]
            subset = np.asarray(subset)

            if ndv is not None:
                subset = subset[subset != ndv * mfac]
            subset = subset[np.isfinite(subset)]
            subset = subset[subset >= -999]
            subset = subset[subset <= 999]
            if subset.size == 0:
                continue

            if mesh.node(i).z() == -2000:
                mean = np.average(subset)
                std = np.std(subset)
                if (mean - 2 * std) < np.min(subset):
                    values[i] = np.min(subset)
                    numvaluesgathered[i] = -2000
                else:
                    subset = subset[subset <= (mean - 2 * std)]
                    values[i] = np.mean(subset)
                    numvaluesgathered[i] = -2000
            else:
                values[i] = values[i] + np.sum(subset)
                numvaluesgathered[i] = numvaluesgathered[i] + np.size(subset)

    return values, numvaluesgathered


def interpolate(mesh, rasterlist, minBathyDepth, mfac, imethod):
    raster_paths = []
    for line in Path(rasterlist).read_text(encoding="utf-8").splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        raster_paths.append(stripped.split()[0])

    val = np.zeros(mesh.numNodes(), dtype=float)
    numval = np.zeros(mesh.numNodes(), dtype=float)
    mesh.size = mesh.computeMeshSize()

    if imethod == "griddata":
        meshconn, xc, yc, boundaryNodes = meshconnectivity(mesh)

    rastersize = 0.0
    N = None
    CA = None
    for raster in raster_paths:
        if imethod == "CA":
            with open_raster(raster) as data:
                newrastersize = get_rastersize(data)
            if abs(rastersize - newrastersize) > 0.10 or N is None or CA is None:
                rastersize = newrastersize
                numCells = np.asarray(compute_numcells(mesh, rastersize))
                N = numCells[0, :]
                CA = numCells[1, :]
            val, numval = gathervalues(mesh, raster, N, CA, mfac, val, numval)
        elif imethod == "griddata":
            val, numval = griddata(mesh, meshconn, xc, yc, boundaryNodes, raster, mfac, val, numval)
        else:
            raise ValueError(f"Unknown interpolation method: {imethod}")

    interpvalues = np.zeros(mesh.numNodes(), dtype=float)
    for i in range(mesh.numNodes()):
        if numval[i] == 0:
            interpvalues[i] = mesh.node(i).z()
        elif numval[i] == -2000:
            interpvalues[i] = val[i]
        else:
            interpvalues[i] = val[i] / numval[i]

        if interpvalues[i] >= 0 and interpvalues[i] < minBathyDepth:
            interpvalues[i] = minBathyDepth

    mesh.setZ(interpvalues)
    return mesh


def compute_numcells(mesh, rastersize):
    sfactor = np.ones(mesh.numNodes(), dtype=float)
    for i in range(mesh.numNodes()):
        if (mesh.node(i).z() < -1001) and (mesh.node(i).z() > -1100):
            sfactor[i] = mesh.node(i).z() * -1 - 1000
        elif mesh.node(i).z() == -2000:
            sfactor[i] = 2.0
        else:
            sfactor[i] = 1.0

    N = (0.25 * np.asarray(mesh.size)) / rastersize
    N = N * sfactor
    N = np.round(N)
    CA = np.where(N < 1, 1, (2 * N + 1) ** 2)
    return N, CA
