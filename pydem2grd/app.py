"""Command-line application for PyDEM2GRD."""

import argparse
from pathlib import Path

from pydem2grd.config import load_configuration
from pydem2grd.mesh import Mesh
from pydem2grd.src.interpolate import interpolate
from pydem2grd.workflow import run_configuration


def _arguments(argv=None):
    parser = argparse.ArgumentParser(
        description="Interpolate DEM elevations onto an ADCIRC fort.14 mesh."
    )
    parser.add_argument(
        "inputs",
        nargs="+",
        metavar="INPUT",
        help="config.json, or INPUT_FORT14 OUTPUT_FORT14 RASTER_LIST",
    )
    parser.add_argument("--multiplication-factor", "-m", type=float, default=-1.0)
    parser.add_argument("--minimum-bathymetric-depth", type=float, default=0.0)
    parser.add_argument("--method", choices=("CA", "griddata"), default="CA")
    return parser.parse_args(argv)


def run(argv=None):
    args = _arguments(argv)
    if len(args.inputs) == 1:
        config_path = Path(args.inputs[0])
        if config_path.suffix.lower() != ".json":
            raise SystemExit("A single input must be a JSON configuration file")
        configuration = load_configuration(config_path)
        print("Reading mesh and validating raster sources...")
        stats, unresolved, report_path = run_configuration(configuration)
        for name, counts in stats.items():
            print(
                "{} nodes flagged: {}\n{} nodes updated: {}\n{} nodes unresolved: {}".format(
                    name,
                    counts["flagged"],
                    name,
                    counts["updated"],
                    name,
                    counts["unresolved"],
                )
            )
        print("Unresolved-node report: {}".format(report_path))
        print("Finished.")
        return
    if len(args.inputs) != 3:
        raise SystemExit("Legacy usage requires INPUT_FORT14 OUTPUT_FORT14 RASTER_LIST")

    input_mesh, output_mesh, raster_list = args.inputs
    mesh = Mesh(input_mesh)
    print("Reading mesh...")
    mesh.read()
    print("Building element table...")
    mesh.buildElementTable()

    print("Interpolating...")
    interpolated_mesh = interpolate(
        mesh,
        raster_list,
        args.minimum_bathymetric_depth,
        args.multiplication_factor,
        args.method,
    )

    print("Saving mesh file...")
    interpolated_mesh.write(output_mesh)

    print("Finished.")
