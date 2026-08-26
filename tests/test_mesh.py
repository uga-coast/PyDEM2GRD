import tempfile
import unittest
from pathlib import Path

from pydem2grd.mesh import Fort14Error, Mesh


FORT14 = """two triangle mesh
2 4
10 0.0 0.0 -1002.0
20 1.0 0.0 -1002.0
30 1.0 1.0 2.5
40 0.0 1.0 -2000.0
7 3 10 20 30
11 3 10 30 40
1 = Number of open boundaries
2 = Total number of open boundary nodes
2 0 = boundary 1
10
20
0 = Number of land boundaries
0 = Total number of land boundary nodes
"""


class MeshTest(unittest.TestCase):
    def _write(self, directory, name="fort.14", contents=FORT14):
        path = Path(directory) / name
        path.write_text(contents, encoding="utf-8")
        return path

    def test_reads_non_contiguous_ids_and_builds_topology(self):
        with tempfile.TemporaryDirectory() as directory:
            mesh = Mesh.from_file(self._write(directory))

        self.assertEqual(mesh.title, "two triangle mesh")
        self.assertEqual(mesh.numNodes(), 4)
        self.assertEqual(mesh.numElements(), 2)
        self.assertEqual(mesh.nodeIndexById(30), 2)
        self.assertEqual(mesh.elementIndexById(11), 1)
        self.assertEqual(mesh.connectivity(), [[10, 20, 30], [10, 30, 40]])
        self.assertEqual(mesh.numElementsAroundNode(mesh.nodeIndexById(10)), 2)
        self.assertEqual({node.id() for node in mesh.boundaryNodes()}, {10, 20, 30, 40})
        for actual, expected in zip(
            mesh.computeMeshSize(), [1.13807119, 1.0, 1.13807119, 1.0]
        ):
            self.assertAlmostEqual(actual, expected)

    def test_round_trip_retains_boundary_section_and_updates_elevations(self):
        with tempfile.TemporaryDirectory() as directory:
            source = self._write(directory)
            output = Path(directory) / "output.14"
            mesh = Mesh.from_file(source)
            mesh.setZ([1.0, 2.0, 3.0, 4.0])
            mesh.write(output)
            result = Mesh.from_file(output)

            self.assertEqual([result.node(i).z() for i in range(4)], [1.0, 2.0, 3.0, 4.0])
            self.assertIn("2 0 = boundary 1", output.read_text(encoding="utf-8"))

    def test_rejects_unknown_element_node(self):
        invalid = FORT14.replace("7 3 10 20 30", "7 3 10 20 99")
        with tempfile.TemporaryDirectory() as directory:
            with self.assertRaisesRegex(Fort14Error, "unknown node ID"):
                Mesh.from_file(self._write(directory, contents=invalid))

    def test_checked_in_example_is_a_valid_mesh(self):
        example = Path(__file__).parents[1] / "example" / "mesh_x1002.grd"
        mesh = Mesh.from_file(example)
        self.assertEqual(mesh.numNodes(), 4596)
        self.assertEqual(mesh.numElements(), 8941)
        self.assertEqual(mesh.node(mesh.nodeIndexById(4599)).id(), 4599)
        self.assertEqual(mesh.elementIndexById(8956), 8940)


if __name__ == "__main__":
    unittest.main()
