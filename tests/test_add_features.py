from RDG import RDG, Node, Edge

import unittest

def test_inserting_orfs_edges():
    g = RDG()
    g = RDG.load_example(g)
    node_key = g.get_new_node_key()

    start_codon_position = 300
    stop_codon_position = 400
    start_node = Node(
        key=node_key,
        node_type="start",
        position=start_codon_position,
        edges_in=[],
        edges_out=[],
        nodes_in=[],
        nodes_out=[],
    )

    node_key = node_key + 1
    stop_node = Node(
        key=node_key,
        node_type="stop",
        position=stop_codon_position,
        edges_in=[],
        edges_out=[],
        nodes_in=[],
        nodes_out=[],
    )
    g.add_node(start_node)
    g.add_node(stop_node)
    g.insert_ORF(g.edges[2], start_node, stop_node)
    assert list(g.edges.keys()) == [1, 2, 3, 4, 5, 6, 7]


def test_insterting_orfs_nodes():
    g = RDG()
    g = RDG.load_example(g)
    node_key = g.get_new_node_key()

    start_codon_position = 15
    stop_codon_position = 25
    start_node = Node(
        key=node_key,
        node_type="start",
        position=start_codon_position,
        edges_in=[],
        edges_out=[],
        nodes_in=[],
        nodes_out=[],
    )

    node_key = node_key + 1
    stop_node = Node(
        key=node_key,
        node_type="stop",
        position=stop_codon_position,
        edges_in=[],
        edges_out=[],
        nodes_in=[],
        nodes_out=[],
    )
    g.add_node(start_node)
    g.add_node(stop_node)
    g.insert_ORF(g.edges[2], start_node, stop_node)
    assert list(g.nodes.keys()) == [1, 2, 3, 4, 5, 6, 7, 8]


def test_insterting_orfs_branch_points():
    g = RDG()
    g = RDG.load_example(g)
    node_key = g.get_new_node_key()

    start_codon_position = 15
    stop_codon_position = 25
    start_node = Node(
        key=node_key,
        node_type="start",
        position=start_codon_position,
        edges_in=[],
        edges_out=[],
        nodes_in=[],
        nodes_out=[],
    )

    node_key = node_key + 1
    stop_node = Node(
        key=node_key,
        node_type="stop",
        position=stop_codon_position,
        edges_in=[],
        edges_out=[],
        nodes_in=[],
        nodes_out=[],
    )
    g.add_node(start_node)
    g.add_node(stop_node)
    g.insert_ORF(g.edges[2], start_node, stop_node)

    branch_points = g.get_branch_points()

    assert len(branch_points) == 2


def invalid_stop_pos():
    g = RDG()
    g = RDG.load_example(g)

    g.add_open_reading_frame(15, 250000)


class TestAddORF(unittest.TestCase):
    def test(self):
        with self.assertRaises(Exception) as context:
            invalid_stop_pos()

        self.assertTrue("Next in frame stop codon" in str(context.exception))


def test_add_open_reading_frame():
    g = RDG()
    g = RDG.load_example(g)

    g.add_open_reading_frame(15, 25)
    branch_points = g.get_branch_points()

    assert len(branch_points) == 2

def invalid_readthrough_stop():
    g = RDG()
    g = RDG.load_example(g)

    g.add_stop_codon_readthrough(15, 250000)

class TestAddSCR(unittest.TestCase):
    def test(self):
        with self.assertRaises(Exception) as context:
            invalid_readthrough_stop()

        self.assertTrue("Next in frame stop codon" in str(context.exception))


def test_add_stop_codon_readthrough():
    g = RDG()
    g = RDG.load_example(g)

    g.add_stop_codon_readthrough(100, 150)
    assert len(g.nodes) == 7

def invalid_frameshift():
    g = RDG()
    g = RDG.load_example(g)

    g.add_frameshift(15, 250000, -1)

class TestAddFS(unittest.TestCase):
    def test(self):
        with self.assertRaises(Exception) as context:
            invalid_frameshift()

        self.assertTrue("Next in frame stop codon" in str(context.exception))


def test_add_frameshift():
    g = RDG()
    g = RDG.load_example(g)

    g.add_frameshift(30, 150, -1)
    assert len(g.nodes) == 8