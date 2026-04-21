from pathlib import Path

from pydem2grd.mesh import read_mesh


def test_read_write_roundtrip(tmp_path: Path):
    mesh_text = """grid
2 4
1 0.0 0.0 -1002.0
2 1.0 0.0 -1002.0
3 1.0 1.0 -1002.0
4 0.0 1.0 -1002.0
1 3 1 2 3
2 3 1 3 4
0 = Number of open boundaries
0 = Total number of open boundary nodes
0 = Number of land boundaries
0 = Total number of land boundary nodes
"""
    inmesh = tmp_path / "mesh.grd"
    outmesh = tmp_path / "mesh_out.grd"
    inmesh.write_text(mesh_text, encoding="utf-8")

    mesh = read_mesh(inmesh)
    assert mesh.numNodes() == 4
    assert mesh.numElements() == 2
    assert mesh.boundary_node_ids() == [1, 2, 3, 4]

    mesh.write(outmesh)
    mesh2 = read_mesh(outmesh)
    assert mesh2.numNodes() == 4
    assert mesh2.numElements() == 2
