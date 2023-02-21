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


def test_add_edge():
    g = RDG()
    g = RDG.load_example(g)

    edge_key = g.get_new_edge_key()
    e = Edge(edge_key, "translated", 2, 1, coordinates=(1, 999))
    g.add_edge(e, 1, 2)
    assert g.nodes[1].output_edges == [1, edge_key]
    assert g.nodes[2].input_edges == [2, edge_key]
