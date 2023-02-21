from RDG import RDG, Node, Edge


def test_edge_init():
    g = RDG()
    edge_key = g.get_new_edge_key()
    edge = Edge(key=edge_key, edge_type="translated", from_node=2, to_node=1, coordinates=(1,999))
    assert edge.key == edge_key
    assert edge.edge_type == "translated"
    assert edge.from_node == 2
    assert edge.to_node == 1
    assert edge.coordinates == (1,999)        

def test_get_frame():
    g = RDG()
    edge_key = g.get_new_edge_key()
    edge = Edge(key=edge_key, edge_type="translated", from_node=2, to_node=1, coordinates=(1,999))
    assert edge.get_frame() == 1
