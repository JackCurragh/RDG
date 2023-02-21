from RDG import RDG, Node, Edge

import unittest

# Set of tests to check the load_example function
nodes = {
    1: Node(key=1, node_type="5_prime", position=0, edges_out=[1], nodes_out=[2]),
    2: Node(key=2, node_type="3_prime", position=999, edges_in=[1], nodes_in=[1]),
}
edges = {1: Edge(1, "untranslated", from_node=1, to_node=2, coordinates=(1, 999))}


def test_load_locus_stop_update():
    g = RDG()
    g = RDG.Load(
        g, locus_name="GeneA", locus_start=0, locus_stop=10000, nodes=nodes, edges=edges
    )
    assert g.locus_stop == 10000


def test_load_nodes():
    g = RDG()
    g = RDG.load_example(g)
    assert len(list(g.nodes.keys())) == 5


def load_graph_no_nodes():
    g = RDG()
    g = RDG.Load(
        g, locus_name="GeneA", locus_start=0, locus_stop=10000, nodes={}, edges=edges
    )


def load_graph_no_edges():
    g = RDG()
    g = RDG.Load(
        g, locus_name="GeneA", locus_start=0, locus_stop=10000, nodes=nodes, edges={}
    )


def stop_greater_than_start_error():
    g = RDG()
    g = RDG.Load(
        g, locus_name="GeneA", locus_start=10, locus_stop=0, nodes=nodes, edges=edges
    )


class TestLoad(unittest.TestCase):
    def test_no_nodes(self):
        with self.assertRaises(Exception) as context:
            load_graph_no_nodes()

        self.assertTrue(
            "nodes and edges must be provided to load a graph" in str(context.exception)
        )

    def test_no_edges(self):
        with self.assertRaises(Exception) as context:
            load_graph_no_edges()

        self.assertTrue(
            "nodes and edges must be provided to load a graph" in str(context.exception)
        )

    def test_stop_greater_than_start(self):
        with self.assertRaises(Exception) as context:
            stop_greater_than_start_error()

        self.assertTrue(
            "locus_stop must be greater than locus_start" in str(context.exception)
        )
