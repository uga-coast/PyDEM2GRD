from __future__ import annotations

import argparse
from pathlib import Path

from pydem2grd.mesh import read_mesh
from pydem2grd.src.interpolate import interpolate


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Interpolate DEM rasters onto an ADCIRC-style unstructured mesh without ADCIRC Modules."
    )
    parser.add_argument("--inmesh", required=True, help="Input ADCIRC .grd / fort.14 mesh")
    parser.add_argument("--outmesh", required=True, help="Output mesh with interpolated z values")
    parser.add_argument("--raster-list", required=True, help="Text file listing raster paths")
    parser.add_argument(
        "--mfac",
        type=float,
        default=-1.0,
        help="Multiplier applied to raster elevations before interpolation. Default: -1.0",
    )
    parser.add_argument(
        "--min-bathy",
        type=float,
        default=0.0,
        help="Minimum non-negative bathymetry depth to enforce after interpolation. Default: 0.0",
    )
    parser.add_argument(
        "--method",
        choices=["CA", "griddata"],
        default="CA",
        help="Interpolation method. Default: CA",
    )
    return parser


def run(argv: list[str] | None = None) -> None:
    args = build_parser().parse_args(argv)

    print(f"Reading mesh: {args.inmesh}")
    mesh = read_mesh(args.inmesh)

    print(f"Interpolating from raster list: {args.raster_list}")
    outmesh = interpolate(mesh, args.raster_list, args.min_bathy, args.mfac, args.method)

    out_path = Path(args.outmesh)
    print(f"Writing mesh: {out_path}")
    outmesh.write(out_path)
    print("Finished.")
