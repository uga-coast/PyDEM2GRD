from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List
import math
import numpy as np


@dataclass(slots=True)
class Node:
    node_id: int
    xcoord: float
    ycoord: float
    elevation: float

    def id(self) -> int:
        return self.node_id

    def x(self) -> float:
        return self.xcoord

    def y(self) -> float:
        return self.ycoord

    def z(self) -> float:
        return self.elevation

    def set_z(self, value: float) -> None:
        self.elevation = float(value)


@dataclass(slots=True)
class Element:
    element_id: int
    node_ids: tuple[int, int, int]

    def id(self) -> int:
        return self.element_id


@dataclass(slots=True)
class Boundary:
    kind: str
    header_tokens: list[str]
    rows: list[list[str]]

    @property
    def node_ids(self) -> list[int]:
        ids: list[int] = []
        for row in self.rows:
            if row:
                try:
                    ids.append(int(row[0]))
                except ValueError:
                    pass
        return ids


class Mesh:
    def __init__(
        self,
        title: str,
        nodes: List[Node],
        elements: List[Element],
        open_boundaries: list[Boundary] | None = None,
        land_boundaries: list[Boundary] | None = None,
    ) -> None:
        self.title = title
        self.nodes = nodes
        self.elements = elements
        self.open_boundaries = open_boundaries or []
        self.land_boundaries = land_boundaries or []
        self._id_to_index = {node.id(): i for i, node in enumerate(nodes)}
        self._node_to_elements: list[list[int]] | None = None
        self._boundary_node_ids: list[int] | None = None
        self.size: np.ndarray | None = None

    @classmethod
    def from_file(cls, path: str | Path) -> "Mesh":
        return read_mesh(path)

    def read(self) -> int:
        return 1

    def write(self, path: str | Path) -> None:
        path = Path(path)
        with path.open("w", encoding="utf-8") as fp:
            fp.write(f"{self.title}\n")
            fp.write(f"{self.numElements()}  {self.numNodes()}\n")
            for node in self.nodes:
                fp.write(
                    f"{node.id():>10d} {node.x():.10f} {node.y():.10f} {node.z():.10f}\n"
                )
            for element in self.elements:
                n1, n2, n3 = element.node_ids
                fp.write(f"{element.id():>10d}    3 {n1:>5d} {n2:>5d} {n3:>5d}\n")

            fp.write(f"{len(self.open_boundaries)} = Number of open boundaries\n")
            fp.write(
                f"{sum(len(b.rows) for b in self.open_boundaries)} = Total number of open boundary nodes\n"
            )
            for boundary in self.open_boundaries:
                header = " ".join(boundary.header_tokens).strip()
                fp.write(f"{header}\n")
                for row in boundary.rows:
                    fp.write(" ".join(row).strip() + "\n")

            fp.write(f"{len(self.land_boundaries)} = Number of land boundaries\n")
            fp.write(
                f"{sum(len(b.rows) for b in self.land_boundaries)} = Total number of land boundary nodes\n"
            )
            for boundary in self.land_boundaries:
                header = " ".join(boundary.header_tokens).strip()
                fp.write(f"{header}\n")
                for row in boundary.rows:
                    fp.write(" ".join(row).strip() + "\n")

    def numNodes(self) -> int:
        return len(self.nodes)

    def numElements(self) -> int:
        return len(self.elements)

    def node(self, index: int) -> Node:
        return self.nodes[index]

    def nodeIndexById(self, node_id: int) -> int:
        return self._id_to_index[node_id]

    def connectivity(self) -> np.ndarray:
        return np.array([element.node_ids for element in self.elements], dtype=int)

    def setZ(self, values: Iterable[float]) -> None:
        for node, value in zip(self.nodes, values):
            node.set_z(float(value))

    def _build_node_to_elements(self) -> list[list[int]]:
        if self._node_to_elements is None:
            mapping: list[list[int]] = [[] for _ in range(self.numNodes())]
            for element_index, element in enumerate(self.elements):
                for node_id in element.node_ids:
                    mapping[self.nodeIndexById(node_id)].append(element_index)
            self._node_to_elements = mapping
        return self._node_to_elements

    def numElementsAroundNode(self, node_index: int) -> int:
        return len(self._build_node_to_elements()[node_index])

    def elementTable(self, node: Node, around_index: int) -> Element:
        node_index = self.nodeIndexById(node.id())
        element_index = self._build_node_to_elements()[node_index][around_index]
        return self.elements[element_index]

    def boundary_node_ids(self) -> list[int]:
        if self._boundary_node_ids is not None:
            return self._boundary_node_ids

        edge_counts: dict[tuple[int, int], int] = {}
        for element in self.elements:
            a, b, c = element.node_ids
            for edge in ((a, b), (b, c), (c, a)):
                key = tuple(sorted(edge))
                edge_counts[key] = edge_counts.get(key, 0) + 1

        node_ids: set[int] = set()
        for (a, b), count in edge_counts.items():
            if count == 1:
                node_ids.add(a)
                node_ids.add(b)

        if not node_ids:
            for boundary in self.open_boundaries + self.land_boundaries:
                node_ids.update(boundary.node_ids)

        self._boundary_node_ids = sorted(node_ids)
        return self._boundary_node_ids

    def boundaryNodes(self) -> list[Node]:
        return [self.nodes[self.nodeIndexById(node_id)] for node_id in self.boundary_node_ids()]

    def computeMeshSize(self) -> np.ndarray:
        node_edges: list[list[float]] = [[] for _ in range(self.numNodes())]
        for element in self.elements:
            ids = element.node_ids
            for a, b in ((ids[0], ids[1]), (ids[1], ids[2]), (ids[2], ids[0])):
                ia = self.nodeIndexById(a)
                ib = self.nodeIndexById(b)
                na = self.node(ia)
                nb = self.node(ib)
                distance = math.hypot(na.x() - nb.x(), na.y() - nb.y())
                node_edges[ia].append(distance)
                node_edges[ib].append(distance)

        sizes = np.zeros(self.numNodes(), dtype=float)
        for i, edges in enumerate(node_edges):
            sizes[i] = float(np.mean(edges)) if edges else 0.0
        self.size = sizes
        return sizes


def _parse_boundary_sections(lines: list[str], start: int, count: int) -> tuple[list[Boundary], int]:
    boundaries: list[Boundary] = []
    idx = start
    for kind_index in range(count):
        header_tokens = lines[idx].split()
        idx += 1
        if not header_tokens:
            raise ValueError(f"Empty boundary header for boundary {kind_index + 1}")
        num_nodes = int(header_tokens[0])
        rows: list[list[str]] = []
        for _ in range(num_nodes):
            rows.append(lines[idx].split())
            idx += 1
        boundaries.append(Boundary(kind="", header_tokens=header_tokens, rows=rows))
    return boundaries, idx


def read_mesh(path: str | Path) -> Mesh:
    path = Path(path)
    with path.open("r", encoding="utf-8", errors="replace") as fp:
        lines = [line.rstrip("\n") for line in fp]

    if len(lines) < 2:
        raise ValueError(f"Mesh file is too short: {path}")

    title = lines[0].strip() or "grid"
    counts = lines[1].split()
    if len(counts) < 2:
        raise ValueError(f"Could not parse element/node counts from: {lines[1]!r}")
    num_elements = int(counts[0])
    num_nodes = int(counts[1])

    idx = 2
    nodes: list[Node] = []
    for _ in range(num_nodes):
        tokens = lines[idx].split()
        idx += 1
        if len(tokens) < 4:
            raise ValueError(f"Invalid node line: {lines[idx - 1]!r}")
        nodes.append(Node(int(tokens[0]), float(tokens[1]), float(tokens[2]), float(tokens[3])))

    elements: list[Element] = []
    for eidx in range(num_elements):
        tokens = lines[idx].split()
        idx += 1
        if len(tokens) < 5:
            raise ValueError(f"Invalid element line: {lines[idx - 1]!r}")
        element_id = int(tokens[0])
        nnodes = int(tokens[1])
        if nnodes != 3:
            raise ValueError(f"Only triangular elements are supported, got {nnodes} on line {idx}")
        elements.append(Element(element_id if element_id > 0 else eidx + 1, (int(tokens[2]), int(tokens[3]), int(tokens[4]))))

    open_boundaries: list[Boundary] = []
    land_boundaries: list[Boundary] = []
    if idx < len(lines):
        open_count = int(lines[idx].split()[0])
        idx += 1
        if idx < len(lines):
            idx += 1  # total open nodes line
        open_boundaries, idx = _parse_boundary_sections(lines, idx, open_count)
        for boundary in open_boundaries:
            boundary.kind = "open"

    if idx < len(lines):
        land_count = int(lines[idx].split()[0])
        idx += 1
        if idx < len(lines):
            idx += 1  # total land nodes line
        land_boundaries, idx = _parse_boundary_sections(lines, idx, land_count)
        for boundary in land_boundaries:
            boundary.kind = "land"

    return Mesh(title=title, nodes=nodes, elements=elements, open_boundaries=open_boundaries, land_boundaries=land_boundaries)
