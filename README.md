# PyDEM2GRD
Interpolate a digital elevation model (DEM) to an ADCIRC unstructured mesh.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

python -m pydem2grd

The following is a list of nodal flag values that are accepted.
* -1000/-1001: Automatic CAA method of Bilskie and Hagen (2012). This flag value will create the most topographically accurate surface.
* -10XX: Flagged values less than -1001 will use a CAA * XX value. This is used for smoothing. For example, a flag value of -1002 will mutliple the default CAA value by 2 thereby increasing the control volume/stencil by 2. A flag value of -1002 is a good choice that balances topographic accuracy and smoothness of the topobathy required for a numerical simulation.
* -2000: This flag is used for vertical/raised feature nodes. Elevation values larger than mean + 2*sigma are averaged so the crown of a feature is captured.

## Docker

There is a Docker image available that contains the necessary prerequisites. You will still need to clone this repository to obtain the latest version of the source code.

https://hub.docker.com/r/mbilskie/pydem2grd/tags

docker pull mbilskie/pydem2grd:1

## Prerequisites

What things you need to install the software and how to install them

* [ADCIRC Modules](https://github.com/zcobell/ADCIRCModules)
* [GDAL](https://pypi.org/project/GDAL/)
* [shapely](https://shapely.readthedocs.io/en/stable/)
* [rasterio](https://rasterio.readthedocs.io/en/stable/#)

### Installing

A step by step series of examples that tell you how to get a development env running

Say what the step will be

```
Give the example
```

And repeat

```
until finished
```

End with an example of getting some data out of the system or using it for a little demo

## Running the tests

Explain how to run the automated tests for this system

### Break down into end to end tests

## Built With

* [ADCIRC Modules](https://github.com/zcobell/ADCIRCModules) - The ADCIRC file I/O used

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

This project is licensed under the *** License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* [Zach Cobell](https://github.com/zcobell)
* Hat tip to anyone whose code was used
* Inspiration
* etc
