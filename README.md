# PyDEM2GRD
Interpolate a digital elevation model (DEM) to an ADCIRC unstructured mesh.

## Getting Started

Create the Conda environment and install the project:

```bash
conda env create -f environment.yml
conda activate pydem2grd
```

Or install it into an existing Python environment:

```bash
python -m pip install -e .
```

Installing the project also provides the equivalent `pydem2grd` command. The
preferred workflow uses a complete JSON configuration:

```bash
pydem2grd config.json
```

The repository includes a runnable example:

```bash
pydem2grd example/config.json
```

Relative paths are resolved from the directory containing `config.json`. A
minimal seamless topo-bathy example is:

```json
{
  "input_mesh": "mesh/fort.14",
  "output_mesh": "output/interpolated.14",
  "crs": "EPSG:26917",
  "elevation_multiplication_factor": 1.0,
  "node_flags": {
    "bathy": -9999,
    "topo": 9999
  },
  "unresolved_report": "output/unresolved.csv",
  "raster_sets": [
    {
      "name": "Seamless topo-bathy DEM",
      "tiles": ["dem/tile_01.tif", "dem/tile_02.tif"],
      "rules": [
        {
          "node_flags": ["bathy"],
          "domain": "bathy",
          "method": "direct_lookup",
          "land_threshold": 0.0,
          "minimum_depth": 1.0,
          "multiplication_factor": -1.0
        },
        {
          "node_flags": ["topo"],
          "domain": "topo",
          "method": "CA",
          "smoothing_factor": 2
        }
      ]
    }
  ]
}
```

Node-flag names are user-defined. Only nodes whose existing z-value equals a
configured sentinel are changed. Ordinary mesh elevations remain untouched,
and unresolved flagged nodes retain their original sentinel.

The original positional interface remains available:

```bash
pydem2grd INPUT_FORT14 OUTPUT_FORT14 RASTER_LIST
```

For example:

```bash
pydem2grd example/mesh_x1002.grd interpolated.grd rasterlist.txt
```

The legacy nodal flags remain supported by that interface:

* `-1000`/`-1001`: automatic CA method.
* `-10XX`: apply the encoded smoothing factor `XX` after establishing the
  integer base radius.
* `-2000`: retain the historical raised-feature processing behavior.

The CLI `--multiplication-factor` converts raster elevation units or sign; it
does not control CA smoothing. Run `pydem2grd --help` for all legacy options.

## Cell-area interpolation

The CA method follows Bilskie and Hagen (2013). Local mesh size
`Delta_M` is the arithmetic mean of the lengths of all unique mesh edges
incident to a node. The continuous base radius is:

```text
N_raw = 0.25 * Delta_M / Delta_DEM
```

If `N_raw < 1`, PyDEM2GRD directly looks up the raster cell containing the
node. Otherwise it converts `N_raw` to an integer using conventional half-up
rounding. The integer base radius is then multiplied by the rule's smoothing
factor. Radius `r` requests a centered `(2r + 1) x (2r + 1)` stencil. A
smoothing factor therefore extends the radius, not the final cell count.

At the outer coverage boundary the stencil remains centered and the
unavailable portion is clipped; it is never shifted or enlarged. Masked,
declared NoData, and nonfinite values are excluded. There is no arbitrary
elevation-range filter.

## Raster priority, tiling, and rules

Raster sets are evaluated in JSON order. Tiles within one set are equal-priority
pieces of one logical DEM and collectively fill a stencil. Overlapping tiles
contribute at most one value per cell; conflicting valid values at the same
aligned cell are reported as an error. Later raster sets fill only locations
left unresolved by earlier sets and never overwrite higher-priority values. A
fallback rule contributes only where its own method and smoothing footprint
overlaps the stencil established by the first contributing raster set.

Every tile in a raster set must use the mesh CRS, square north-up pixels, one
resolution, and an aligned grid. Lower-priority sets may fill individual gaps
only when their grids are compatible with the stencil's established grid. An
incompatible lower-resolution or shifted grid is rejected if it would be
needed to complete a partial stencil. If a higher-priority set contributes no
values at all, the next eligible set may establish a new stencil using its own
native resolution. Automatic reprojection and resampling are intentionally not
performed.

Rules associate node-flag classes with per-raster-set behavior. `topo` rules
use every otherwise valid elevation, including values below datum. `bathy`
rules accept only raster values below `land_threshold` (default `0.0`). Minimum
depth is applied only to bathymetric results. A rule may override the global
elevation multiplication factor.

At completion, PyDEM2GRD reports flagged, updated, and unresolved counts for
each node-flag class. It also writes an unresolved-node CSV beside the output
mesh unless `unresolved_report` specifies another path. The CSV includes node
ID, coordinates, flag class, original sentinel, and reason.

## Docker

There is a Docker image available that contains the necessary prerequisites. You will still need to clone this repository to obtain the latest version of the source code.

https://hub.docker.com/r/mbilskie/pydem2grd/tags

docker pull mbilskie/pydem2grd:1

## Dependencies

* [NumPy](https://numpy.org/)
* [Shapely](https://shapely.readthedocs.io/en/stable/)
* [Rasterio](https://rasterio.readthedocs.io/en/stable/)

ADCIRC `fort.14`/`.grd` files are read and written natively. Boundary sections
are retained when a mesh is written, while boundary nodes used by interpolation
are derived from the element topology.

```python
from pydem2grd import Mesh

mesh = Mesh.from_file("fort.14")
mesh.write("fort-updated.14")
```

## Running the tests

```bash
python -m unittest discover -s tests
```

## Built With

* Native Python ADCIRC `fort.14` mesh I/O
* Rasterio-based DEM access (no direct `osgeo.gdal` dependency)

## Contributing

* In progress...

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Authors

* **Matthew V Bilskie, PhD** - *Initial work* - (https://github.com/mattbilskie)

See also the list of [contributors](https://github.com/mattbilskie/PyDEM2GRD/contributors) who participated in this project.

## Citation

Please appropriatley cite this work in publications, techincal reports, source code, etc. as:

```
@article{Bilskie:2015,
   author = {Bilskie, Matthew V. and Coggin, David and Hagen, Scott C. and Medeiros, Stephen C.},
   title = {Terrain-driven unstructured mesh development through semi-automatic vertical feature extraction},
   journal = {Advances in Water Resources},
   volume = {86, Part A},
   pages = {102-118},
   ISSN = {0309-1708},
   DOI = {http://dx.doi.org/10.1016/j.advwatres.2015.09.020},
   url = {http://www.sciencedirect.com/science/article/pii/S0309170815002274},
   year = {2015},
   type = {Journal Article}
}
```

## Notes

This project is still under development.

## License

This project is licensed under the GNU General Public License v3.0; see
[`LICENSE`](LICENSE) for details.

## Acknowledgments

* [Zach Cobell](https://github.com/zcobell)
* Hat tip to anyone whose code was used
* Inspiration
* etc
