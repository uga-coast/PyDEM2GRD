"""Native ADCIRC ``fort.14`` mesh support.

The small compatibility API in this module is intentionally shaped like the
mesh operations that PyDEM2GRD uses. This keeps the interpolation code focused
on interpolation while making mesh I/O a local, testable part of the project.
"""

from collections import Counter
from dataclasses import dataclass
from math import hypot
from pathlib import Path


class Fort14Error(ValueError):
    """Raised when an ADCIRC mesh is malformed or unsupported."""


@dataclass
class Node:
    _id: int
    _x: float
    _y: float
    _z: float

    def id(self):
        return self._id

    def x(self):
        return self._x

    def y(self):
        return self._y

    def z(self):
        return self._z


@dataclass
class Element:
    _id: int
    node_ids: tuple

    def id(self):
        return self._id


class _ElementTableBuilder:
    def __init__(self, mesh):
        self._mesh = mesh

    def build(self):
        self._mesh.buildElementTable()


class _Topology:
    def __init__(self, mesh):
        self._mesh = mesh

    def elementTable(self):
        return _ElementTableBuilder(self._mesh)


class Mesh:
    """An in-memory ADCIRC mesh with native ``fort.14`` read/write support."""

    def __init__(self, filename=None):
        self.filename = Path(filename) if filename is not None else None
        self.title = ""
        self._nodes = []
        self._elements = []
        self._node_indices = {}
        self._element_indices = {}
        self._element_table = None
        self._boundary_lines = []
        self.size = []

    @classmethod
    def from_file(cls, filename):
        mesh = cls(filename)
        mesh.read()
        return mesh

    def read(self, filename=None):
        if filename is not None:
            self.filename = Path(filename)
        if self.filename is None:
            raise ValueError("A fort.14 filename is required")

        try:
            with self.filename.open("r", encoding="utf-8") as stream:
                lines = stream.readlines()
        except OSError as exc:
            raise Fort14Error("Could not read {}: {}".format(self.filename, exc)) from exc

        if len(lines) < 2:
            raise Fort14Error("{} must contain a title and size header".format(self.filename))

        self.title = lines[0].rstrip("\r\n")
        counts = lines[1].split()
        if len(counts) < 2:
            raise Fort14Error("Line 2 must contain element and node counts")
        try:
            element_count, node_count = int(counts[0]), int(counts[1])
        except ValueError as exc:
            raise Fort14Error("Line 2 contains invalid element or node counts") from exc
        if element_count < 0 or node_count < 0:
            raise Fort14Error("Element and node counts cannot be negative")

        cursor = 2
        nodes = []
        for _ in range(node_count):
            line_number = cursor + 1
            if cursor >= len(lines):
                raise Fort14Error("Unexpected end of file while reading nodes")
            fields = lines[cursor].split()
            cursor += 1
            if len(fields) < 4:
                raise Fort14Error("Line {} is not a valid node record".format(line_number))
            try:
                nodes.append(Node(int(fields[0]), float(fields[1]), float(fields[2]), float(fields[3])))
            except ValueError as exc:
                raise Fort14Error("Line {} contains an invalid node record".format(line_number)) from exc

        node_ids = {node.id() for node in nodes}
        if len(node_ids) != len(nodes):
            raise Fort14Error("Node IDs must be unique")

        elements = []
        for _ in range(element_count):
            line_number = cursor + 1
            if cursor >= len(lines):
                raise Fort14Error("Unexpected end of file while reading elements")
            fields = lines[cursor].split()
            cursor += 1
            try:
                element_id = int(fields[0])
                vertex_count = int(fields[1])
                element_nodes = tuple(int(value) for value in fields[2:])
            except (IndexError, ValueError) as exc:
                raise Fort14Error("Line {} is not a valid element record".format(line_number)) from exc
            if vertex_count < 3 or len(element_nodes) != vertex_count:
                raise Fort14Error(
                    "Line {} declares {} element nodes but contains {}".format(
                        line_number, vertex_count, len(element_nodes)
                    )
                )
            missing = set(element_nodes) - node_ids
            if missing:
                raise Fort14Error(
                    "Line {} references unknown node ID(s): {}".format(
                        line_number, ", ".join(str(value) for value in sorted(missing))
                    )
                )
            elements.append(Element(element_id, element_nodes))

        if len({element.id() for element in elements}) != len(elements):
            raise Fort14Error("Element IDs must be unique")

        self._nodes = nodes
        self._elements = elements
        self._node_indices = {node.id(): index for index, node in enumerate(nodes)}
        self._element_indices = {element.id(): index for index, element in enumerate(elements)}
        self._boundary_lines = [line.rstrip("\r\n") for line in lines[cursor:]]
        self._element_table = None
        self.size = []
        return True

    def write(self, filename=None):
        destination = Path(filename) if filename is not None else self.filename
        if destination is None:
            raise ValueError("A fort.14 output filename is required")

        with destination.open("w", encoding="utf-8", newline="\n") as stream:
            stream.write("{}\n".format(self.title))
            stream.write("{} {}\n".format(self.numElements(), self.numNodes()))
            for node in self._nodes:
                stream.write(
                    "{} {:.16g} {:.16g} {:.16g}\n".format(
                        node.id(), node.x(), node.y(), node.z()
                    )
                )
            for element in self._elements:
                stream.write(
                    "{} {} {}\n".format(
                        element.id(), len(element.node_ids), " ".join(map(str, element.node_ids))
                    )
                )
            boundary_lines = self._boundary_lines or [
                "0 = Number of open boundaries",
                "0 = Total number of open boundary nodes",
                "0 = Number of land boundaries",
                "0 = Total number of land boundary nodes",
            ]
            for line in boundary_lines:
                stream.write("{}\n".format(line))

        self.filename = destination

    def numNodes(self):
        return len(self._nodes)

    def numElements(self):
        return len(self._elements)

    def node(self, index):
        return self._nodes[index]

    def nodeIndexById(self, node_id):
        try:
            return self._node_indices[int(node_id)]
        except KeyError as exc:
            raise KeyError("Unknown node ID {}".format(node_id)) from exc

    def elementIndexById(self, element_id):
        try:
            return self._element_indices[int(element_id)]
        except KeyError as exc:
            raise KeyError("Unknown element ID {}".format(element_id)) from exc

    def connectivity(self):
        return [list(element.node_ids) for element in self._elements]

    def buildElementTable(self):
        table = [[] for _ in self._nodes]
        for element in self._elements:
            for node_id in element.node_ids:
                table[self.nodeIndexById(node_id)].append(element)
        self._element_table = table

    def topology(self):
        return _Topology(self)

    def elementTable(self, node, index):
        if self._element_table is None:
            self.buildElementTable()
        return self._element_table[self.nodeIndexById(node.id())][index]

    def numElementsAroundNode(self, node_index):
        if self._element_table is None:
            self.buildElementTable()
        return len(self._element_table[node_index])

    def boundaryNodes(self):
        edge_counts = Counter()
        for element in self._elements:
            node_ids = element.node_ids
            for index, node_id in enumerate(node_ids):
                edge_counts[tuple(sorted((node_id, node_ids[(index + 1) % len(node_ids)])))] += 1
        boundary_ids = {
            node_id
            for edge, count in edge_counts.items()
            if count == 1
            for node_id in edge
        }
        return [node for node in self._nodes if node.id() in boundary_ids]

    def computeMeshSize(self):
        lengths = [[] for _ in self._nodes]
        seen_edges = set()
        for element in self._elements:
            node_ids = element.node_ids
            for index, first_id in enumerate(node_ids):
                second_id = node_ids[(index + 1) % len(node_ids)]
                edge = tuple(sorted((first_id, second_id)))
                if edge in seen_edges:
                    continue
                seen_edges.add(edge)
                first = self.node(self.nodeIndexById(first_id))
                second = self.node(self.nodeIndexById(second_id))
                length = hypot(second.x() - first.x(), second.y() - first.y())
                lengths[self.nodeIndexById(first_id)].append(length)
                lengths[self.nodeIndexById(second_id)].append(length)
        return [sum(values) / len(values) if values else 0.0 for values in lengths]

    def setZ(self, values):
        if len(values) != self.numNodes():
            raise ValueError("Expected {} elevations, got {}".format(self.numNodes(), len(values)))
        for node, value in zip(self._nodes, values):
            node._z = float(value)
