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

```bash
python -m pydem2grd INPUT_FORT14 OUTPUT_FORT14 RASTER_LIST
```

Installing the project also provides the equivalent `pydem2grd` command.

For example:

```bash
python -m pydem2grd example/mesh_x1002.grd interpolated.grd rasterlist.txt
```

Run `python -m pydem2grd --help` for interpolation method, multiplication
factor, and minimum-depth options.

The following is a list of nodal flag values that are accepted.

* -1000/-1001: Automatic CA method of Bilskie and Hagen (2012). This flag value will create the most topographically accurate surface.
* -10XX: Flagged values less than -1001 use a CA scale factor of XX. This is used for smoothing. For example, -1002 multiplies the default control-area radius by 2, increasing the interpolation stencil.
* -2000: This flag is used for vertical/raised feature nodes. Elevation values larger than mean + 2*sigma are averaged so the crown of a feature is captured.

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
