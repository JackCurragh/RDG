

from RDG import RDG, Node, Edge

import unittest

def test_newick():
    g = RDG()
    g = RDG.load_example(g)
    newick = g.newick()
    assert newick == "((2:990,(5:900)4:90)3:10)1;"

def test_newick_spare_endpoint():
    g = RDG()
    g = RDG.load_example(g)
    g.add_open_reading_frame(30, 40)

    newick = g.newick()

    assert newick == "(((5:900)4:90,((8:960)7:10,2:970)6:20)3:10)1;"

def newick_error_too_many_roots():
    g = RDG()
    g = RDG.load_example(g)
    g.add_node(Node(6, "5_prime", 10, [], [], [], []))
    newick = g.newick()


class TestLoad(unittest.TestCase):
    def test_newick(self):
        with self.assertRaises(Exception) as context:
            newick_error_too_many_roots()

        self.assertTrue(
            "Graph has more than one startpoint" in str(context.exception)
        )

