import json
import tempfile
import unittest
from pathlib import Path

import numpy as np
import rasterio
from rasterio.transform import from_origin

from pydem2grd.config import load_configuration
from pydem2grd.mesh import Mesh
from pydem2grd.workflow import ca_radius, run_configuration, round_half_up


class ConfiguredWorkflowTest(unittest.TestCase):
    crs = "EPSG:26917"

    def _raster(self, path, values, origin_x=0.0, origin_y=3.0, resolution=1.0,
                nodata=-9999.0, crs=None, y_resolution=None):
        y_resolution = resolution if y_resolution is None else y_resolution
        with rasterio.open(
            path,
            "w",
            driver="GTiff",
            height=values.shape[0],
            width=values.shape[1],
            count=1,
            dtype=values.dtype,
            transform=from_origin(origin_x, origin_y, resolution, y_resolution),
            nodata=nodata,
            crs=crs or self.crs,
        ) as dataset:
            dataset.write(values, 1)
        return path

    def _single_node_mesh(self, path, x=1.5, y=1.5, elevation=9999.0):
        path.write_text(
            "single node\n0 1\n1 {} {} {}\n".format(x, y, elevation)
            + "0 = Number of open boundaries\n"
            + "0 = Total number of open boundary nodes\n"
            + "0 = Number of land boundaries\n"
            + "0 = Total number of land boundary nodes\n",
            encoding="utf-8",
        )
        return path

    def _radius_one_mesh(self, path, elevation=9999.0):
        path.write_text(
            "radius one\n3 4\n"
            "1 1.5 1.5 {}\n".format(elevation)
            + "2 5.5 1.5 0\n3 1.5 5.5 0\n4 -2.5 1.5 0\n"
            + "1 3 1 2 3\n2 3 1 3 4\n3 3 1 4 2\n"
            + "0 = Number of open boundaries\n"
            + "0 = Total number of open boundary nodes\n"
            + "0 = Number of land boundaries\n"
            + "0 = Total number of land boundary nodes\n",
            encoding="utf-8",
        )
        return path

    def _config(self, directory, mesh, raster_sets, flags=None, **extra):
        config = {
            "input_mesh": Path(mesh).name,
            "output_mesh": "output.14",
            "crs": self.crs,
            "node_flags": flags or {"topo": 9999},
            "raster_sets": raster_sets,
        }
        config.update(extra)
        path = Path(directory) / "config.json"
        path.write_text(json.dumps(config), encoding="utf-8")
        return load_configuration(path)

    def _set(self, tiles, rules=None, name="DEM"):
        return {
            "name": name,
            "tiles": [Path(tile).name for tile in tiles],
            "rules": rules or [
                {"node_flags": ["topo"], "domain": "topo", "method": "CA"}
            ],
        }

    def test_canonical_ca_radii(self):
        self.assertEqual([ca_radius(size, 5.0) for size in (20, 40, 80, 160)], [1, 2, 4, 8])

    def test_raw_radius_below_one_is_always_direct_lookup(self):
        self.assertEqual(ca_radius(19.999, 5.0, smoothing_factor=10), 0)

    def test_half_up_rounding_and_smoothing_order(self):
        self.assertEqual(round_half_up(2.5), 3)
        self.assertEqual(round_half_up(1.5), 2)
        self.assertEqual(ca_radius(50.0, 5.0), 3)
        self.assertEqual(ca_radius(32.0, 5.0, smoothing_factor=2), 4)

    def test_relative_paths_and_direct_lookup(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            mesh = self._single_node_mesh(directory / "mesh.14")
            raster = self._raster(
                directory / "dem.tif", np.arange(1, 10, dtype="float32").reshape(3, 3)
            )
            config = self._config(
                directory,
                mesh,
                [self._set([raster], [{"node_flags": ["topo"], "domain": "topo", "method": "direct_lookup"}])],
            )
            stats, unresolved, report = run_configuration(config)
            self.assertEqual(Mesh.from_file(directory / "output.14").node(0).z(), 5.0)
            self.assertEqual(stats["topo"], {"flagged": 1, "updated": 1, "unresolved": 0})
            self.assertEqual(unresolved, [])
            self.assertTrue(report.exists())

    def test_adjacent_tiles_complete_centered_stencil(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            mesh = self._radius_one_mesh(directory / "mesh.14")
            values = np.arange(1, 13, dtype="float32").reshape(3, 4)
            left = self._raster(directory / "left.tif", values[:, :2], origin_x=0)
            right = self._raster(directory / "right.tif", values[:, 2:], origin_x=2)
            config = self._config(directory, mesh, [self._set([left, right])])
            run_configuration(config)
            self.assertEqual(Mesh.from_file(directory / "output.14").node(0).z(), 6.0)

    def test_stencil_is_clipped_not_shifted_at_coverage_edge(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            mesh = self._radius_one_mesh(directory / "mesh.14")
            raster = self._raster(
                directory / "dem.tif", np.arange(1, 10, dtype="float32").reshape(3, 3), origin_x=1
            )
            config = self._config(directory, mesh, [self._set([raster])])
            run_configuration(config)
            # Only stencil columns centered at x=1.5 and x=2.5 are available.
            self.assertEqual(Mesh.from_file(directory / "output.14").node(0).z(), 4.5)

    def test_nodata_falls_back_without_overwriting_priority_values(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            mesh = self._radius_one_mesh(directory / "mesh.14")
            high_values = np.full((3, 3), 10, dtype="float32")
            high_values[1, 1] = -9999
            high = self._raster(directory / "high.tif", high_values)
            low = self._raster(directory / "low.tif", np.full((3, 3), 20, dtype="float32"))
            config = self._config(
                directory, mesh,
                [self._set([high], name="high"), self._set([low], name="low")],
            )
            run_configuration(config)
            self.assertAlmostEqual(
                Mesh.from_file(directory / "output.14").node(0).z(), 100.0 / 9.0
            )

    def test_lower_direct_lookup_rule_does_not_fill_outer_ca_cells(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            mesh = self._radius_one_mesh(directory / "mesh.14")
            high_values = np.full((3, 3), -9999, dtype="float32")
            high_values[1, 1] = 10
            high = self._raster(directory / "high.tif", high_values)
            low = self._raster(directory / "low.tif", np.full((3, 3), 20, dtype="float32"))
            ca_rule = [{"node_flags": ["topo"], "domain": "topo", "method": "CA"}]
            direct_rule = [
                {"node_flags": ["topo"], "domain": "topo", "method": "direct_lookup"}
            ]
            config = self._config(
                directory,
                mesh,
                [self._set([high], ca_rule, "high"), self._set([low], direct_rule, "low")],
            )
            run_configuration(config)
            self.assertEqual(Mesh.from_file(directory / "output.14").node(0).z(), 10.0)

    def test_overlapping_tiles_are_not_double_weighted(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            mesh = self._radius_one_mesh(directory / "mesh.14")
            first = self._raster(
                directory / "first.tif", np.arange(1, 10, dtype="float32").reshape(3, 3)
            )
            second = self._raster(
                directory / "second.tif", np.arange(1, 10, dtype="float32").reshape(3, 3)
            )
            config = self._config(directory, mesh, [self._set([first, second])])
            run_configuration(config)
            self.assertEqual(Mesh.from_file(directory / "output.14").node(0).z(), 5.0)

    def test_valid_same_set_tile_fills_overlapping_nodata(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            mesh = self._single_node_mesh(directory / "mesh.14")
            nodata = self._raster(
                directory / "nodata.tif", np.full((3, 3), -9999, dtype="float32")
            )
            valid = self._raster(
                directory / "valid.tif", np.full((3, 3), 7, dtype="float32")
            )
            rules = [{"node_flags": ["topo"], "domain": "topo", "method": "direct_lookup"}]
            config = self._config(directory, mesh, [self._set([nodata, valid], rules)])
            run_configuration(config)
            self.assertEqual(Mesh.from_file(directory / "output.14").node(0).z(), 7.0)

    def test_unresolved_flag_is_preserved_and_unflagged_node_is_untouched(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            mesh = directory / "mesh.14"
            mesh.write_text(
                "unresolved\n0 2\n1 20 20 9999\n2 1.5 1.5 123.5\n0\n0\n0\n0\n",
                encoding="utf-8",
            )
            raster = self._raster(directory / "dem.tif", np.ones((3, 3), dtype="float32"))
            config = self._config(directory, mesh, [self._set([raster])])
            stats, unresolved, report = run_configuration(config)
            output = Mesh.from_file(directory / "output.14")
            self.assertEqual(output.node(0).z(), 9999.0)
            self.assertEqual(output.node(1).z(), 123.5)
            self.assertEqual(stats["topo"]["unresolved"], 1)
            self.assertEqual(unresolved[0]["node_id"], 1)
            self.assertIn("no valid raster value", report.read_text(encoding="utf-8"))

    def test_seamless_dem_applies_asymmetric_bathy_and_topo_rules(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            mesh_path = directory / "mesh.14"
            mesh_path.write_text(
                "flags\n0 2\n1 0.5 1.5 -9999\n2 1.5 0.5 9999\n"
                "0\n0\n0\n0\n", encoding="utf-8"
            )
            raster = self._raster(
                directory / "dem.tif", np.array([[0, 0], [-2, -4]], dtype="float32"),
                origin_y=2,
            )
            rules = [
                {"node_flags": ["bathy"], "domain": "bathy", "method": "direct_lookup", "minimum_depth": 3, "multiplication_factor": -1},
                {"node_flags": ["topo"], "domain": "topo", "method": "direct_lookup"},
            ]
            config = self._config(
                directory, mesh_path, [self._set([raster], rules)],
                flags={"bathy": -9999, "topo": 9999},
            )
            run_configuration(config)
            output = Mesh.from_file(directory / "output.14")
            self.assertEqual(output.node(0).z(), -9999.0)  # land-side zero is ineligible
            self.assertEqual(output.node(1).z(), -4.0)     # topo keeps water-side data

    def test_minimum_depth_applies_only_to_bathy(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            mesh_path = directory / "mesh.14"
            mesh_path.write_text(
                "flags\n0 2\n1 0.5 0.5 -9999\n2 1.5 0.5 9999\n0\n0\n0\n0\n",
                encoding="utf-8",
            )
            raster = self._raster(
                directory / "dem.tif", np.array([[-0.25, 0.25]], dtype="float32"), origin_y=1
            )
            rules = [
                {"node_flags": ["bathy"], "domain": "bathy", "method": "direct_lookup", "minimum_depth": 1, "multiplication_factor": -1},
                {"node_flags": ["topo"], "domain": "topo", "method": "direct_lookup", "minimum_depth": 1},
            ]
            config = self._config(directory, mesh_path, [self._set([raster], rules)], flags={"bathy": -9999, "topo": 9999})
            run_configuration(config)
            output = Mesh.from_file(directory / "output.14")
            self.assertEqual(output.node(0).z(), 1.0)
            self.assertEqual(output.node(1).z(), 0.25)

    def test_crs_mismatch_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            mesh = self._single_node_mesh(directory / "mesh.14")
            raster = self._raster(directory / "dem.tif", np.ones((3, 3), dtype="float32"), crs="EPSG:4326")
            config = self._config(directory, mesh, [self._set([raster])])
            with self.assertRaisesRegex(ValueError, "CRS mismatch"):
                run_configuration(config)

    def test_non_square_pixels_are_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            mesh = self._single_node_mesh(directory / "mesh.14")
            raster = self._raster(
                directory / "dem.tif", np.ones((3, 3), dtype="float32"), y_resolution=2
            )
            config = self._config(directory, mesh, [self._set([raster])])
            with self.assertRaisesRegex(ValueError, "non-square"):
                run_configuration(config)

    def test_incompatible_grid_cannot_fill_partial_stencil(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            mesh = self._radius_one_mesh(directory / "mesh.14")
            high_values = np.ones((3, 3), dtype="float32")
            high_values[1, 1] = -9999
            high = self._raster(directory / "high.tif", high_values)
            low = self._raster(
                directory / "low.tif", np.ones((6, 6), dtype="float32"),
                origin_y=3, resolution=0.5,
            )
            config = self._config(directory, mesh, [self._set([high]), self._set([low])])
            with self.assertRaisesRegex(ValueError, "incompatible grid"):
                run_configuration(config)

    def test_incompatible_lower_grid_can_establish_stencil_when_higher_set_is_empty(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            mesh = self._single_node_mesh(directory / "mesh.14")
            high = self._raster(
                directory / "high.tif", np.full((3, 3), -9999, dtype="float32")
            )
            low = self._raster(
                directory / "low.tif", np.full((6, 6), 8, dtype="float32"),
                origin_y=3, resolution=0.5,
            )
            rules = [{"node_flags": ["topo"], "domain": "topo", "method": "direct_lookup"}]
            config = self._config(
                directory,
                mesh,
                [self._set([high], rules, "high"), self._set([low], rules, "low")],
            )
            run_configuration(config)
            self.assertEqual(Mesh.from_file(directory / "output.14").node(0).z(), 8.0)


if __name__ == "__main__":
    unittest.main()
