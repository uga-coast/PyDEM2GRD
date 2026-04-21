# PyDEM2GRD

A clean, installable rewrite of PyDEM2GRD that removes the ADCIRC Modules dependency while preserving the ADCIRC-style mesh workflow.

## Features

- Pure-Python ADCIRC `.grd` / `fort.14` mesh reader and writer
- No `pyadcircmodules` dependency
- Console command: `pydem2grd`
- Supports `CA` and `griddata` interpolation methods
- Uses `rasterio` and `shapely` for raster and geometry handling

## Installation

Using pip in an existing environment:

```bash
pip install -e .
```

Then verify:

```bash
pydem2grd --help
```

## Usage

```bash
pydem2grd \
  --inmesh mesh.grd \
  --outmesh mesh_interp.grd \
  --raster-list rasterlist.txt \
  --mfac -1.0 \
  --min-bathy 0.0 \
  --method CA
```

## Required inputs

- `--inmesh`: input ADCIRC `.grd` or `fort.14` mesh
- `--outmesh`: output mesh path
- `--raster-list`: text file with one raster path per line

## Notes

- The input mesh workflow is designed around flagged Z values, similar to the original code.
- `CA` remains the safest default.
- Mesh and rasters should use the same coordinate reference system.
