"""Configured, priority-aware DEM interpolation workflow."""

import csv
from dataclasses import dataclass
from decimal import Decimal, ROUND_HALF_UP
from math import floor, isclose

import numpy as np
import rasterio
from rasterio.crs import CRS
from rasterio.windows import Window

from .mesh import Mesh


def round_half_up(value):
    """Return conventional nearest-integer, half-up rounding for a nonnegative value."""
    if value < 0:
        raise ValueError("CA radius cannot be negative")
    return int(Decimal(str(float(value))).quantize(Decimal("1"), rounding=ROUND_HALF_UP))


def ca_radius(mesh_size, raster_size, smoothing_factor=1):
    """Compute the final integer CA radius using the Bilskie-Hagen relationship."""
    if raster_size <= 0:
        raise ValueError("Raster size must be positive")
    if smoothing_factor < 1:
        raise ValueError("Smoothing factor must be positive")
    raw_radius = 0.25 * float(mesh_size) / float(raster_size)
    if raw_radius < 1:
        return 0
    return round_half_up(raw_radius) * int(smoothing_factor)


@dataclass(frozen=True)
class Grid:
    crs: CRS
    resolution: float
    origin_x: float
    origin_y: float

    def compatible(self, other):
        tolerance = max(self.resolution, other.resolution) * 1e-8
        if self.crs != other.crs or not isclose(
            self.resolution, other.resolution, rel_tol=1e-9, abs_tol=tolerance
        ):
            return False
        x_offset = (other.origin_x - self.origin_x) / self.resolution
        y_offset = (other.origin_y - self.origin_y) / self.resolution
        return isclose(x_offset, round(x_offset), abs_tol=1e-7) and isclose(
            y_offset, round(y_offset), abs_tol=1e-7
        )

    def containing_cell(self, x, y):
        return (
            int(floor((self.origin_y - y) / self.resolution)),
            int(floor((x - self.origin_x) / self.resolution)),
        )

    def center(self, row, col):
        return (
            self.origin_x + (col + 0.5) * self.resolution,
            self.origin_y - (row + 0.5) * self.resolution,
        )


class Tile:
    def __init__(self, path, mesh_crs):
        self.path = path
        self.dataset = rasterio.open(path)
        try:
            self.grid = _validated_grid(self.dataset, mesh_crs, path)
        except Exception:
            self.dataset.close()
            raise

    def close(self):
        self.dataset.close()

    def values_at(self, points, rule, multiplication_factor):
        """Read requested aligned cell centers with one bounded raster window."""
        if not points:
            return {}
        keys = list(points)
        coordinates = np.asarray([points[key] for key in keys], dtype=float)
        x = coordinates[:, 0]
        y = coordinates[:, 1]
        rows = np.floor(
            (self.grid.origin_y - y) / self.grid.resolution + 1e-9
        ).astype(int)
        cols = np.floor(
            (x - self.grid.origin_x) / self.grid.resolution + 1e-9
        ).astype(int)
        inside = (
            (rows >= 0)
            & (cols >= 0)
            & (rows < self.dataset.height)
            & (cols < self.dataset.width)
        )
        if not np.any(inside):
            return {}
        selected = np.flatnonzero(inside)
        selected_rows = rows[selected]
        selected_cols = cols[selected]
        row_start = int(selected_rows.min())
        col_start = int(selected_cols.min())
        row_stop = int(selected_rows.max()) + 1
        col_stop = int(selected_cols.max()) + 1
        window = Window(col_start, row_start, col_stop - col_start, row_stop - row_start)
        raster_values = self.dataset.read(1, window=window, masked=True).astype(float)
        values = raster_values[selected_rows - row_start, selected_cols - col_start]
        valid = ~np.ma.getmaskarray(values) & np.isfinite(np.asarray(values))
        raw_values = np.asarray(values)
        if rule.domain == "bathy":
            valid &= raw_values < rule.land_threshold
        factor = (
            rule.multiplication_factor
            if rule.multiplication_factor is not None
            else multiplication_factor
        )
        return {
            keys[selected[position]]: float(raw_values[position]) * factor
            for position in np.flatnonzero(valid)
        }


class RasterProduct:
    def __init__(self, raster_set, mesh_crs):
        self.definition = raster_set
        self.tiles = []
        try:
            for path in raster_set.tiles:
                self.tiles.append(Tile(path, mesh_crs))
            self.grid = self.tiles[0].grid
            for tile in self.tiles[1:]:
                if not self.grid.compatible(tile.grid):
                    raise ValueError(
                        "Raster set {!r} contains tiles on incompatible grids: {}".format(
                            raster_set.name, tile.path
                        )
                    )
        except Exception:
            self.close()
            raise

    def close(self):
        for tile in self.tiles:
            tile.close()

    def rule_for(self, flag_name):
        for rule in self.definition.rules:
            if flag_name in rule.node_flags:
                return rule
        return None

    def values_at(self, points, rule, factor):
        values = {}
        for tile in self.tiles:
            tile_values = tile.values_at(points, rule, factor)
            for key, value in tile_values.items():
                if key in values and not isclose(
                    values[key], value, rel_tol=1e-9, abs_tol=1e-9
                ):
                    x, y = points[key]
                    raise ValueError(
                        "Raster set {!r} has conflicting valid values in overlapping "
                        "tiles at ({}, {})".format(self.definition.name, x, y)
                    )
                values[key] = value
        return values


def _validated_grid(dataset, mesh_crs, path):
    if dataset.crs is None:
        raise ValueError("Raster {} does not declare a CRS".format(path))
    if dataset.crs != mesh_crs:
        raise ValueError(
            "CRS mismatch: mesh uses {}, but raster {} uses {}".format(
                mesh_crs, path, dataset.crs
            )
        )
    transform = dataset.transform
    tolerance = max(abs(transform.a), abs(transform.e), 1.0) * 1e-10
    if abs(transform.b) > tolerance or abs(transform.d) > tolerance:
        raise ValueError("Raster {} has rotated or sheared pixels".format(path))
    if transform.a <= 0 or transform.e >= 0:
        raise ValueError("Raster {} is not a north-up raster grid".format(path))
    x_resolution = abs(float(transform.a))
    y_resolution = abs(float(transform.e))
    if not isclose(x_resolution, y_resolution, rel_tol=1e-9, abs_tol=tolerance):
        raise ValueError(
            "Raster {} has non-square pixels ({} x {})".format(
                path, x_resolution, y_resolution
            )
        )
    if x_resolution <= 0:
        raise ValueError("Raster {} has an invalid pixel size".format(path))
    return Grid(dataset.crs, x_resolution, float(transform.c), float(transform.f))


def _stencil(grid, x, y, radius):
    center_row, center_col = grid.containing_cell(x, y)
    return {
        (row_offset, col_offset): grid.center(
            center_row + row_offset, center_col + col_offset
        )
        for row_offset in range(-radius, radius + 1)
        for col_offset in range(-radius, radius + 1)
    }


def _rule_radius(rule, mesh_size, grid):
    if rule.method == "direct_lookup":
        return 0
    return ca_radius(mesh_size, grid.resolution, rule.smoothing_factor)


def _candidate_values(product, target, rule, factor, unresolved_only=None):
    if unresolved_only is None:
        points = target
    else:
        points = {key: target[key] for key in unresolved_only}
    return product.values_at(points, rule, factor)


def _apply_minimum_depth(value, rule, configuration):
    if rule.domain != "bathy":
        return value
    minimum = (
        rule.minimum_depth
        if rule.minimum_depth is not None
        else configuration.minimum_bathymetric_depth
    )
    if minimum > 0 and 0 <= value < minimum:
        return minimum
    return value


def interpolate_configured(mesh, configuration):
    """Interpolate configured node flags using ordered, priority-mosaiced raster sets."""
    try:
        mesh_crs = CRS.from_user_input(configuration.crs)
    except Exception as exc:
        raise ValueError("Invalid mesh CRS {!r}: {}".format(configuration.crs, exc)) from exc

    products = []
    try:
        for item in configuration.raster_sets:
            products.append(RasterProduct(item, mesh_crs))
        mesh_sizes = mesh.computeMeshSize()
        results = [node.z() for node in (mesh.node(i) for i in range(mesh.numNodes()))]
        stats = {
            name: {"flagged": 0, "updated": 0, "unresolved": 0}
            for name in configuration.node_flags
        }
        unresolved = []

        sentinel_to_name = {value: name for name, value in configuration.node_flags.items()}
        for index in range(mesh.numNodes()):
            node = mesh.node(index)
            flag_name = sentinel_to_name.get(node.z())
            if flag_name is None:
                continue
            stats[flag_name]["flagged"] += 1
            target = None
            target_grid = None
            target_rule = None
            gathered = {}

            for product in products:
                rule = product.rule_for(flag_name)
                if rule is None:
                    continue
                if target is None:
                    radius = _rule_radius(rule, mesh_sizes[index], product.grid)
                    candidate_target = _stencil(product.grid, node.x(), node.y(), radius)
                    candidate = _candidate_values(
                        product, candidate_target, rule, configuration.multiplication_factor
                    )
                    if not candidate:
                        continue
                    target = candidate_target
                    target_grid = product.grid
                    target_rule = rule
                    gathered.update(candidate)
                else:
                    missing = target.keys() - gathered.keys()
                    if not missing:
                        break
                    product_radius = _rule_radius(rule, mesh_sizes[index], product.grid)
                    product_footprint = _stencil(
                        product.grid, node.x(), node.y(), product_radius
                    )
                    eligible_missing = missing.intersection(product_footprint)
                    if not eligible_missing:
                        continue
                    candidate = _candidate_values(
                        product,
                        target,
                        rule,
                        configuration.multiplication_factor,
                        eligible_missing,
                    )
                    if not candidate:
                        continue
                    if not target_grid.compatible(product.grid):
                        raise ValueError(
                            "Raster set {!r} would fill a partial stencil for node {} using "
                            "an incompatible grid".format(product.definition.name, node.id())
                        )
                    gathered.update(candidate)
                if len(gathered) == len(target):
                    break

            if gathered:
                value = float(np.mean(list(gathered.values())))
                results[index] = _apply_minimum_depth(value, target_rule, configuration)
                stats[flag_name]["updated"] += 1
            else:
                stats[flag_name]["unresolved"] += 1
                unresolved.append(
                    {
                        "node_id": node.id(),
                        "x": node.x(),
                        "y": node.y(),
                        "node_flag": flag_name,
                        "sentinel": node.z(),
                        "reason": "no valid raster value",
                    }
                )

        mesh.setZ(results)
        return stats, unresolved
    finally:
        for product in products:
            product.close()


def write_unresolved_report(path, rows):
    """Write a machine-readable unresolved-node CSV report."""
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = ("node_id", "x", "y", "node_flag", "sentinel", "reason")
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def run_configuration(configuration):
    """Read, interpolate, report, and write one configured mesh run."""
    mesh = Mesh.from_file(configuration.input_mesh)
    stats, unresolved = interpolate_configured(mesh, configuration)
    configuration.output_mesh.parent.mkdir(parents=True, exist_ok=True)
    mesh.write(configuration.output_mesh)
    report_path = configuration.unresolved_report or configuration.output_mesh.with_name(
        configuration.output_mesh.name + ".unresolved.csv"
    )
    write_unresolved_report(report_path, unresolved)
    return stats, unresolved, report_path
