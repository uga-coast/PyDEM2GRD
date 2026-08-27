import tempfile
import unittest
from pathlib import Path

try:
    import numpy as np
    import rasterio
    from rasterio.transform import from_origin

    from pydem2grd.mesh import Mesh
    from pydem2grd.src.interpolate import interpolate
    from pydem2grd.src.raster import (
        coord2pixel,
        get_boundingbox,
        get_numrowcol,
        get_rastersize,
        pixel2coord,
    )
except ImportError:
    np = None


@unittest.skipIf(np is None, "scientific dependencies are not installed")
class InterpolationTest(unittest.TestCase):
    def _raster(self, directory, values, nodata=-9999.0):
        path = Path(directory) / "dem.tif"
        with rasterio.open(
            path,
            "w",
            driver="GTiff",
            height=values.shape[0],
            width=values.shape[1],
            count=1,
            dtype=values.dtype,
            transform=from_origin(0.0, 3.0, 1.0, 1.0),
            nodata=nodata,
        ) as dataset:
            dataset.write(values, 1)
        return path

    def _mesh(self, directory, elevation=-1001.0):
        path = Path(directory) / "mesh.14"
        path.write_text(
            "single node\n"
            "0 1\n"
            "10 1.5 1.5 {}\n".format(elevation)
            + "0 = Number of open boundaries\n"
            + "0 = Total number of open boundary nodes\n"
            + "0 = Number of land boundaries\n"
            + "0 = Total number of land boundary nodes\n",
            encoding="utf-8",
        )
        return Mesh.from_file(path)

    def _square_mesh(self, directory):
        path = Path(directory) / "square.14"
        path.write_text(
            "square mesh\n"
            "2 4\n"
            "10 0.5 0.5 -1001\n"
            "20 2.5 0.5 -1001\n"
            "30 2.5 2.5 -1001\n"
            "40 0.5 2.5 -1001\n"
            "7 3 10 20 30\n"
            "11 3 10 30 40\n"
            "0 = Number of open boundaries\n"
            "0 = Total number of open boundary nodes\n"
            "0 = Number of land boundaries\n"
            "0 = Total number of land boundary nodes\n",
            encoding="utf-8",
        )
        return Mesh.from_file(path)

    def test_raster_coordinate_helpers(self):
        with tempfile.TemporaryDirectory() as directory:
            raster_path = self._raster(
                directory, np.arange(1, 10, dtype="float32").reshape(3, 3)
            )
            with rasterio.open(raster_path) as dataset:
                self.assertEqual(get_boundingbox(dataset), [0.0, 0.0, 3.0, 3.0])
                self.assertEqual(get_numrowcol(dataset), (3, 3))
                self.assertEqual(get_rastersize(dataset), 1.0)
                self.assertEqual(coord2pixel(1.5, 1.5, dataset), (1, 1))
                self.assertEqual(pixel2coord(1, 1, dataset), (1.0, 2.0))

    def test_ca_interpolation_reads_center_cell_for_zero_radius(self):
        with tempfile.TemporaryDirectory() as directory:
            raster_path = self._raster(
                directory, np.arange(1, 10, dtype="float32").reshape(3, 3)
            )
            raster_list = Path(directory) / "rasters.txt"
            raster_list.write_text(str(raster_path) + "\n", encoding="utf-8")

            mesh = interpolate(self._mesh(directory), raster_list, 0.0, 1.0, "CA")
            self.assertEqual(mesh.node(0).z(), 5.0)

    def test_ca_interpolation_ignores_nodata(self):
        with tempfile.TemporaryDirectory() as directory:
            values = np.arange(1, 10, dtype="float32").reshape(3, 3)
            values[1, 1] = -9999.0
            raster_path = self._raster(directory, values)
            raster_list = Path(directory) / "rasters.txt"
            raster_list.write_text(str(raster_path) + "\n", encoding="utf-8")

            mesh = interpolate(self._mesh(directory), raster_list, 0.0, 1.0, "CA")
            self.assertEqual(mesh.node(0).z(), -1001.0)

    def test_ca_interpolation_does_not_reject_finite_values_by_range(self):
        with tempfile.TemporaryDirectory() as directory:
            values = np.full((3, 3), 1500.0, dtype="float32")
            raster_path = self._raster(directory, values)
            raster_list = Path(directory) / "rasters.txt"
            raster_list.write_text(str(raster_path) + "\n", encoding="utf-8")

            mesh = interpolate(self._mesh(directory), raster_list, 0.0, 1.0, "CA")
            self.assertEqual(mesh.node(0).z(), 1500.0)

    def test_griddata_supports_non_contiguous_element_ids(self):
        with tempfile.TemporaryDirectory() as directory:
            raster_path = self._raster(
                directory, np.arange(1, 10, dtype="float32").reshape(3, 3)
            )
            raster_list = Path(directory) / "rasters.txt"
            raster_list.write_text(str(raster_path) + "\n", encoding="utf-8")

            mesh = interpolate(
                self._square_mesh(directory), raster_list, 0.0, 1.0, "griddata"
            )
            self.assertTrue(
                any(mesh.node(index).z() != -1001.0 for index in range(mesh.numNodes()))
            )

    def test_checked_in_example_ca_smoke(self):
        repository = Path(__file__).parents[1]
        mesh = Mesh.from_file(repository / "example" / "mesh_x1002.grd")
        original = [mesh.node(index).z() for index in range(mesh.numNodes())]

        with tempfile.TemporaryDirectory() as directory:
            raster_list = Path(directory) / "rasters.txt"
            raster_list.write_text(
                str(repository / "example" / "Example_DEM_GA.tif") + "\n",
                encoding="utf-8",
            )
            interpolate(mesh, raster_list, 0.0, -1.0, "CA")

        changed = sum(
            mesh.node(index).z() != original[index] for index in range(mesh.numNodes())
        )
        self.assertGreater(changed, 0)


if __name__ == "__main__":
    unittest.main()
