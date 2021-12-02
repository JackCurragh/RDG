from RDG import RDG, Node, Edge


def test_basic_edges():
    assert len(RDG().edges) == 1

def test_basic_nodes():
    assert list(RDG().nodes.keys()) == [1, 2]