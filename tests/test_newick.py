

from RDG import RDG, Node, Edge

import unittest

def test_newick():
    g = RDG()
    g = RDG.load_example(g)
    newick = g.newick()
    assert newick == "((2:1,5:1)3:1);"

def test_newick_spare_endpoint():
    g = RDG()
    g = RDG.load_example(g)
    g.add_open_reading_frame(30, 40)

    newick = g.newick()

    assert newick == "(((2:1,8:1)6:1,5:2)3:2);"

def newick_error_too_many_roots():
    g = RDG()
    g = RDG.load_example(g)
    g.add_node(Node(6, "5_prime", 10, [], [], [], []))
    newick = g.newick()

def newick_error_mutlibranch():
    g = RDG()
    g = RDG.load_example(g)
    g.add_node(Node(6, "3_prime", 10, [], [], [], []))
    g.add_edge(Edge(10, "translated", 3, 6, coordinates=(10,20)), 3, 6)
    newick = g.newick()

class TestLoad(unittest.TestCase):
    def test_newick(self):
        with self.assertRaises(Exception) as context:
            newick_error_too_many_roots()

        self.assertTrue(
            "Graph has more than one startpoint" in str(context.exception)
        )

    def test_newick_multibranch(self):
        with self.assertRaises(Exception) as context:
            newick_error_mutlibranch()

        self.assertTrue(
            "Branchpoint has more than two endpoints" in str(context.exception)
        )