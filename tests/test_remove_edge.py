from RDG import RDG, Node, Edge


def test_remove_edge_from_node():
    g = RDG()
    g = RDG.load_example(g)
    g.remove_edge(2)
    assert g.nodes[3].output_edges == [3]

def test_remove_edge_to_node():
    g = RDG()
    g = RDG.load_example(g)
    g.remove_edge(2)
    assert g.nodes[2].input_edges == []

def test_remove_edge_edges():
    g = RDG()
    g = RDG.load_example(g)
    g.remove_edge(2)
    tf = 2 in g.edges
    assert tf == False