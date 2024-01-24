from RDG import RDG, Node, Edge

import unittest


def test_get_nodes():
    g = RDG()
    g = RDG.load_example(g)
    assert g.get_node_keys() == [1, 2, 3, 4, 5]


def test_get_edges():
    g = RDG()
    g = RDG.load_example(g)
    assert g.get_edge_keys() == [1, 2, 3, 4]


def test_get_edges_from_to():
    g = RDG()
    g = RDG.load_example(g)
    edges = g.get_edges_from_to()
    assert sorted(edges) == sorted([(1, 3), (3, 4), (4, 5), (3, 2)])

def test_remove_node():
    g = RDG()
    g = RDG.load_example(g)
    g.remove_node(3)
    assert g.get_node_keys() == [1, 2, 4, 5]

def test_get_new_node_key():
    g = RDG()
    g = RDG.load_example(g)
    assert g.get_new_node_key() == 6


def test_get_new_node_key_empty_graph():
    g = RDG()
    g.remove_node(1)
    g.remove_node(2)
    assert g.get_new_node_key() == 1


def test_get_new_edge_key():
    g = RDG()
    g = RDG.load_example(g)
    assert g.get_new_edge_key() == 5

def test_get_new_edge_key_empty_graph():
    g = RDG()
    g.remove_node(1)
    g.remove_node(2)
    assert g.get_new_edge_key() == 1

def test_get_key_from_position():
    g = RDG()
    g = RDG.load_example(g)
    assert g.get_key_from_position(10, node_type="start")[0] == 3

def get_key_from_position_error_no_type():
    g = RDG()
    g = RDG.load_example(g)
    g.get_key_from_position(10, node_type="readthrough")

def get_key_from_position_error_no_pos():
    g = RDG()
    g = RDG.load_example(g)
    g.get_key_from_position(10000, node_type="stop")

class TestLoad(unittest.TestCase):
    def test_no_nodes_of_x_type(self):
        with self.assertRaises(Exception) as context:
            get_key_from_position_error_no_type()

        self.assertTrue(
            "There are no nodes of type" in str(context.exception)
        )

    def test_no_nodes_at_pos(self):
        with self.assertRaises(Exception) as context:
            get_key_from_position_error_no_pos()
        
        
        self.assertTrue(
            "There is no node of type 'stop' at position 10000" in str(context.exception)
        )

def test_get_startpoints():
    g = RDG()
    g = RDG.load_example(g)
    startpoints = g.get_startpoints()
    assert startpoints == [1]

def test_get_endpoints():
    g = RDG()
    g = RDG.load_example(g)
    endpoints = g.get_endpoints()
    assert endpoints == [2, 5]

def test_get_start_nodes():
    g = RDG()
    g = RDG.load_example(g)
    start_nodes = g.get_start_nodes()
    assert start_nodes == [3]

def test_get_stop_nodes():
    g = RDG()
    g = RDG.load_example(g)
    stop_nodes = g.get_stop_nodes()
    assert stop_nodes == [4]

def test_is_input_edge_translated():
    g = RDG()
    g = RDG.load_example(g)
    assert g.is_input_edge_translated(4) == True
    assert g.is_input_edge_translated(3) == False

def test_get_orfs_simple():
    g = RDG()
    g = RDG.load_example(g)
    orfs = g.get_orfs()
    assert orfs == [(10, 100)]

def test_get_orfs_w_frameshift():
    g = RDG()
    g = RDG.load_example(g)
    g.add_frameshift(30, 40, 1)
    orfs = g.get_orfs()
    assert orfs == [(10, 30, 100), (10, 30, 40)]

def test_get_orfs_scr():
    g = RDG()
    g = RDG.load_example(g)
    g.add_stop_codon_readthrough(100, 110)
    orfs = g.get_orfs()
    assert orfs == [(10, 100), (10, 110)]

def test_get_frameshifts():
    g = RDG()
    g = RDG.load_example(g)
    g.add_frameshift(30, 40, 1)
    frameshifts = g.get_frameshifts()
    assert frameshifts == [6]

def test_get_unique_paths():
    g = RDG()
    g = RDG.load_example(g)
    paths = g.get_unique_paths()
    assert paths == [[1,3,2], [1,3,4,5]]

def test_get_upstream_branchpoint():
    g = RDG()
    g = RDG.load_example(g)
    bp = g.get_upstream_branchpoint(3)
    assert bp == 1

def test_get_upstream_branchpoint_no_direct_bp():
    g = RDG()
    g = RDG.load_example(g)
    bp = g.get_upstream_branchpoint(5)
    assert bp == 3