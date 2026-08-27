"""Command-line application for PyDEM2GRD."""

import argparse

from pydem2grd.mesh import Mesh
from pydem2grd.src.interpolate import interpolate


def _arguments(argv=None):
    parser = argparse.ArgumentParser(
        description="Interpolate DEM elevations onto an ADCIRC fort.14 mesh."
    )
    parser.add_argument("input_mesh", help="input fort.14/ADCIRC grid file")
    parser.add_argument("output_mesh", help="output fort.14/ADCIRC grid file")
    parser.add_argument("raster_list", help="text file containing one raster path per line")
    parser.add_argument("--multiplication-factor", "-m", type=float, default=-1.0)
    parser.add_argument("--minimum-bathymetric-depth", type=float, default=0.0)
    parser.add_argument("--method", choices=("CA", "griddata"), default="CA")
    return parser.parse_args(argv)


def run(argv=None):
    args = _arguments(argv)
    mesh = Mesh(args.input_mesh)
    print("Reading mesh...")
    mesh.read()
    print("Building element table...")
    mesh.buildElementTable()

    print("Interpolating...")
    interpolated_mesh = interpolate(
        mesh,
        args.raster_list,
        args.minimum_bathymetric_depth,
        args.multiplication_factor,
        args.method,
    )
    
    print("Saving mesh file...")
    interpolated_mesh.write(args.output_mesh)

    print("Finished.")
