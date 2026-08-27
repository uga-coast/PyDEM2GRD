"""DEM-to-mesh interpolation routines."""

import math
import operator
from functools import reduce
from pathlib import Path

import numpy as np
import rasterio
from rasterio.mask import mask
from shapely.geometry import Point, Polygon, box, mapping

from .raster import (
    coord2pixel,
    get_boundingbox,
    get_numrowcol,
    get_rastersize,
    open_raster,
    read_band_as_array,
)


def griddata(
    mesh,
    meshconn,
    xc,
    yc,
    boundary_nodes,
    raster,
    multiplication_factor,
    values,
    num_values_gathered,
):
    """Gather raster cells within a centroid-based control polygon."""
    mesh.size = mesh.computeMeshSize()

    with open_raster(raster) as data:
        bbox = get_boundingbox(data)
        bbox_polygon = box(*bbox)
        raster_size = get_rastersize(data)

    for index in range(mesh.numNodes()):
        node = mesh.node(index)
        buffer_radius = 1.25 * mesh.size[index] * raster_size
        if not bbox_polygon.buffer(buffer_radius).contains(Point(node.x(), node.y())):
            continue
        if node.z() > -999.0:
            continue

        num_elements = mesh.numElementsAroundNode(index)
        point_list = []
        for surrounding_index in range(num_elements):
            element = mesh.elementTable(node, surrounding_index)
            element_index = mesh.elementIndexById(element.id())
            point_list.append((xc[element_index], yc[element_index]))

        if node.id() in boundary_nodes:
            if not bbox_polygon.contains(Point(node.x(), node.y())):
                continue
            point_list.append((node.x(), node.y()))

            if num_elements == 1:
                element = mesh.elementTable(node, 0)
                element_index = mesh.elementIndexById(element.id())
                for neighbor_id in meshconn[element_index]:
                    if neighbor_id == node.id():
                        continue
                    neighbor = mesh.node(mesh.nodeIndexById(neighbor_id))
                    point_list.append(
                        (
                            0.5 * (node.x() + neighbor.x()),
                            0.5 * (node.y() + neighbor.y()),
                        )
                    )

        if len(point_list) < 3:
            continue

        center = tuple(
            map(
                operator.truediv,
                reduce(lambda first, second: list(map(operator.add, first, second)), point_list),
                [len(point_list)] * 2,
            )
        )
        point_list = sorted(
            point_list,
            key=lambda coord: (
                -135
                - math.degrees(
                    math.atan2(*tuple(map(operator.sub, coord, center))[::-1])
                )
            )
            % 360,
        )
        voronoi_polygon = Polygon(point_list)
        if not voronoi_polygon.is_valid or voronoi_polygon.area == 0:
            continue

        with rasterio.open(raster) as source:
            raster_polygon = box(*source.bounds)
            if not raster_polygon.intersects(voronoi_polygon):
                continue
            out_image, _ = mask(
                source,
                [mapping(voronoi_polygon)],
                crop=True,
                filled=False,
            )

        subset = out_image.compressed() * multiplication_factor
        subset = subset[np.isfinite(subset)]
        if subset.size == 0:
            continue
        num_values_gathered[index] += subset.size
        values[index] += np.sum(subset)

    return values, num_values_gathered


def meshconnectivity(mesh):
    """Return element connectivity, centroids, and topological boundary IDs."""
    meshconn = mesh.connectivity()
    xc = np.zeros(mesh.numElements())
    yc = np.zeros(mesh.numElements())

    for index, node_ids in enumerate(meshconn):
        if len(node_ids) != 3:
            raise ValueError("griddata interpolation requires triangular elements")
        nodes = [mesh.node(mesh.nodeIndexById(node_id)) for node_id in node_ids]
        xc[index] = sum(node.x() for node in nodes) / 3.0
        yc[index] = sum(node.y() for node in nodes) / 3.0

    boundary_nodes = [int(node.id()) for node in mesh.boundaryNodes()]
    return meshconn, xc, yc, boundary_nodes


def gathervalues(
    mesh,
    raster,
    cell_radii,
    control_areas,
    multiplication_factor,
    values,
    num_values_gathered,
):
    """Gather raster values from square control areas around mesh nodes."""
    with open_raster(raster) as data:
        num_cols, num_rows = get_numrowcol(data)
        raster_values = read_band_as_array(data, masked=True).astype(float)
        raster_values *= multiplication_factor
        raster_size = get_rastersize(data)
        bbox = get_boundingbox(data)

        for index in range(mesh.numNodes()):
            node = mesh.node(index)
            if node.z() > -999.0:
                continue

            radius = int(cell_radii[index])
            coordinate_radius = (radius + 1) * raster_size
            if (
                node.x() < bbox[0] - coordinate_radius
                or node.x() > bbox[2] + coordinate_radius
                or node.y() < bbox[1] - coordinate_radius
                or node.y() > bbox[3] + coordinate_radius
            ):
                continue

            if num_values_gathered[index] == control_areas[index]:
                continue

            col, row = coord2pixel(node.x(), node.y(), data)
            if col < 0 or row < 0:
                continue

            left = max(col - radius, 0)
            right = min(col + radius + 1, num_cols)
            top = max(row - radius, 0)
            bottom = min(row + radius + 1, num_rows)
            if left >= right or top >= bottom:
                continue

            subset = raster_values[top:bottom, left:right].compressed()
            subset = subset[np.isfinite(subset)]
            subset = subset[(subset >= -999) & (subset <= 999)]
            if subset.size == 0:
                continue

            if node.z() == -2000:
                mean = np.average(subset)
                standard_deviation = np.std(subset)
                threshold = mean - 2 * standard_deviation
                if threshold < np.min(subset):
                    values[index] = np.min(subset)
                else:
                    raised_values = subset[subset <= threshold]
                    values[index] = (
                        np.mean(raised_values) if raised_values.size else np.min(subset)
                    )
                num_values_gathered[index] = -2000
            else:
                values[index] += np.sum(subset)
                num_values_gathered[index] += subset.size

    return values, num_values_gathered


def _raster_paths(raster_list):
    raster_list = Path(raster_list)
    base_directory = raster_list.resolve().parent
    paths = []
    for line in raster_list.read_text(encoding="utf-8").splitlines():
        value = line.strip()
        if value and not value.startswith("#"):
            path = Path(value.split()[0])
            paths.append(str(path if path.is_absolute() else base_directory / path))
    return paths


def interpolate(mesh, raster_list, minimum_bathymetric_depth, multiplication_factor, method):
    """Interpolate every flagged mesh node from the listed DEM rasters."""
    if method not in ("CA", "griddata"):
        raise ValueError("Unknown interpolation method: {}".format(method))

    raster_paths = _raster_paths(raster_list)
    if not raster_paths:
        raise ValueError("Raster list does not contain any raster paths")

    values = np.zeros(mesh.numNodes())
    num_values = np.zeros(mesh.numNodes())
    mesh.size = mesh.computeMeshSize()

    if method == "griddata":
        print("Computing mesh connectivity...")
        meshconn, xc, yc, boundary_nodes = meshconnectivity(mesh)
        print("done!")

    raster_size = None
    cell_radii = None
    control_areas = None
    for raster in raster_paths:
        print(raster)
        if method == "CA":
            with open_raster(raster) as data:
                new_raster_size = get_rastersize(data)
            if raster_size is None or abs(raster_size - new_raster_size) > 0.10:
                print(
                    "Raster size changed from {} to {}. Re-calculating control areas.".format(
                        raster_size or 0.0, new_raster_size
                    )
                )
                raster_size = new_raster_size
                cell_radii, control_areas = compute_numcells(mesh, raster_size)

            values, num_values = gathervalues(
                mesh,
                raster,
                cell_radii,
                control_areas,
                multiplication_factor,
                values,
                num_values,
            )
        else:
            values, num_values = griddata(
                mesh,
                meshconn,
                xc,
                yc,
                boundary_nodes,
                raster,
                multiplication_factor,
                values,
                num_values,
            )

    interpolated_values = np.zeros(mesh.numNodes())
    for index in range(mesh.numNodes()):
        if num_values[index] == 0:
            interpolated_values[index] = mesh.node(index).z()
        elif num_values[index] == -2000:
            interpolated_values[index] = values[index]
        else:
            interpolated_values[index] = values[index] / num_values[index]

        if 0 <= interpolated_values[index] < minimum_bathymetric_depth:
            interpolated_values[index] = minimum_bathymetric_depth

    mesh.setZ(interpolated_values)
    return mesh


def compute_numcells(mesh, raster_size):
    """Compute raster-cell radii and expected control-area sizes per node."""
    if raster_size <= 0:
        raise ValueError("Raster size must be positive")

    scale_factors = np.ones(mesh.numNodes())
    for index in range(mesh.numNodes()):
        elevation = mesh.node(index).z()
        if -1100 < elevation < -1001:
            scale_factors[index] = -elevation - 1000
        elif elevation == -2000:
            scale_factors[index] = 2.0

    cell_radii = (0.25 * np.asarray(mesh.size)) / raster_size
    cell_radii = np.round(cell_radii * scale_factors)
    control_areas = np.where(cell_radii < 1, 1, (2 * cell_radii + 1) ** 2)
    return cell_radii, control_areas
