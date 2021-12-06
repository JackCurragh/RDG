from RDG import RDG, Node, Edge


def test_basic_edges():
    assert len(RDG().edges) == 1

def test_basic_nodes():
    assert list(RDG().nodes.keys()) == [1, 2]

def test_root_to_node_acyclic_node_path():
    g = RDG()
    g = RDG.load_example(g)
    assert g.root_to_node_of_acyclic_node_path(5) == [1, 3, 4, 5]