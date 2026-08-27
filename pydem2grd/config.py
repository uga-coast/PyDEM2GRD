"""JSON configuration loading for the full PyDEM2GRD workflow."""

from dataclasses import dataclass
import json
from pathlib import Path


@dataclass(frozen=True)
class RasterRule:
    node_flags: tuple
    domain: str
    method: str
    smoothing_factor: int = 1
    land_threshold: float = 0.0
    minimum_depth: float = None
    multiplication_factor: float = None


@dataclass(frozen=True)
class RasterSet:
    name: str
    tiles: tuple
    rules: tuple


@dataclass(frozen=True)
class Configuration:
    source: Path
    input_mesh: Path
    output_mesh: Path
    crs: str
    node_flags: dict
    raster_sets: tuple
    multiplication_factor: float = 1.0
    minimum_bathymetric_depth: float = 0.0
    unresolved_report: Path = None


def _path(value, base, field):
    if not isinstance(value, str) or not value.strip():
        raise ValueError("{} must be a non-empty path".format(field))
    path = Path(value)
    return path if path.is_absolute() else base / path


def _number(value, field):
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError("{} must be a number".format(field))
    return float(value)


def _parse_flags(raw):
    if not isinstance(raw, dict) or not raw:
        raise ValueError("node_flags must be a non-empty object")
    flags = {}
    used = {}
    for name, raw_value in raw.items():
        if not isinstance(name, str) or not name.strip():
            raise ValueError("node flag names must be non-empty strings")
        value = raw_value.get("value") if isinstance(raw_value, dict) else raw_value
        sentinel = _number(value, "node_flags.{}".format(name))
        if sentinel in used:
            raise ValueError(
                "node flags {} and {} use the same sentinel {}".format(
                    used[sentinel], name, sentinel
                )
            )
        flags[name] = sentinel
        used[sentinel] = name
    return flags


def _parse_rule(raw, flag_names, set_name, index):
    if not isinstance(raw, dict):
        raise ValueError("raster set {!r} rule {} must be an object".format(set_name, index))
    names = raw.get("node_flags")
    if not isinstance(names, list) or not names:
        raise ValueError("raster set {!r} rule {} requires node_flags".format(set_name, index))
    unknown = set(names) - set(flag_names)
    if unknown:
        raise ValueError("unknown node flag(s): {}".format(", ".join(sorted(unknown))))
    if len(set(names)) != len(names):
        raise ValueError("node_flags in a rule must not contain duplicates")

    domain = str(raw.get("domain", "topo")).lower()
    if domain not in ("topo", "bathy"):
        raise ValueError("rule domain must be 'topo' or 'bathy'")
    method = str(raw.get("method", "CA")).lower()
    if method not in ("ca", "direct_lookup"):
        raise ValueError("rule method must be 'CA' or 'direct_lookup'")
    smoothing = raw.get("smoothing_factor", 1)
    if isinstance(smoothing, bool) or not isinstance(smoothing, int) or smoothing < 1:
        raise ValueError("smoothing_factor must be a positive integer")
    raw_minimum_depth = raw.get("minimum_depth")
    minimum_depth = (
        _number(raw_minimum_depth, "minimum_depth")
        if raw_minimum_depth is not None
        else None
    )
    if minimum_depth is not None and minimum_depth < 0:
        raise ValueError("minimum_depth cannot be negative")
    multiplication = raw.get("multiplication_factor")
    if multiplication is not None:
        multiplication = _number(multiplication, "multiplication_factor")
    return RasterRule(
        tuple(names),
        domain,
        method,
        smoothing,
        _number(raw.get("land_threshold", 0.0), "land_threshold"),
        minimum_depth,
        multiplication,
    )


def load_configuration(filename):
    """Load, validate, and resolve paths in a JSON configuration."""
    source = Path(filename).resolve()
    try:
        raw = json.loads(source.read_text(encoding="utf-8"))
    except OSError as exc:
        raise ValueError("Could not read configuration {}: {}".format(source, exc)) from exc
    except json.JSONDecodeError as exc:
        raise ValueError("Invalid JSON configuration {}: {}".format(source, exc)) from exc
    if not isinstance(raw, dict):
        raise ValueError("Configuration root must be a JSON object")

    base = source.parent
    mesh = raw.get("mesh", {})
    if not isinstance(mesh, dict):
        raise ValueError("mesh must be an object")
    input_mesh = raw.get("input_mesh", mesh.get("input"))
    output_mesh = raw.get("output_mesh", mesh.get("output"))
    crs = raw.get("crs", mesh.get("crs"))
    if not isinstance(crs, str) or not crs.strip():
        raise ValueError("crs must be a non-empty CRS string")

    flags = _parse_flags(raw.get("node_flags"))
    raw_sets = raw.get("raster_sets")
    if not isinstance(raw_sets, list) or not raw_sets:
        raise ValueError("raster_sets must be a non-empty array")
    raster_sets = []
    for set_index, raw_set in enumerate(raw_sets):
        if not isinstance(raw_set, dict):
            raise ValueError("raster_sets[{}] must be an object".format(set_index))
        name = raw_set.get("name", "Raster set {}".format(set_index + 1))
        if not isinstance(name, str) or not name.strip():
            raise ValueError("raster set name must be a non-empty string")
        raw_tiles = raw_set.get("tiles")
        if not isinstance(raw_tiles, list) or not raw_tiles:
            raise ValueError("raster set {!r} requires at least one tile".format(name))
        tiles = tuple(_path(value, base, "raster tile") for value in raw_tiles)
        raw_rules = raw_set.get("rules")
        if not isinstance(raw_rules, list) or not raw_rules:
            raise ValueError("raster set {!r} requires at least one rule".format(name))
        rules = tuple(
            _parse_rule(rule, flags, name, index)
            for index, rule in enumerate(raw_rules)
        )
        seen = set()
        for rule in rules:
            overlap = seen.intersection(rule.node_flags)
            if overlap:
                raise ValueError(
                    "raster set {!r} has multiple rules for {}".format(
                        name, ", ".join(sorted(overlap))
                    )
                )
            seen.update(rule.node_flags)
        raster_sets.append(RasterSet(name, tiles, rules))

    unresolved = raw.get("unresolved_report")
    minimum_depth = _number(
        raw.get("minimum_bathymetric_depth", 0.0), "minimum_bathymetric_depth"
    )
    if minimum_depth < 0:
        raise ValueError("minimum_bathymetric_depth cannot be negative")
    return Configuration(
        source=source,
        input_mesh=_path(input_mesh, base, "input_mesh"),
        output_mesh=_path(output_mesh, base, "output_mesh"),
        crs=crs,
        node_flags=flags,
        raster_sets=tuple(raster_sets),
        multiplication_factor=_number(
            raw.get("elevation_multiplication_factor", 1.0),
            "elevation_multiplication_factor",
        ),
        minimum_bathymetric_depth=minimum_depth,
        unresolved_report=(
            _path(unresolved, base, "unresolved_report") if unresolved else None
        ),
    )
