from RDG import RDG, Node, Edge

def test_get_vertices():
    g = RDG()
    g = RDG.load_example(g)
    assert g.vertices() == [1,2,3,4,5]

def test_get_edges():
    g = RDG()
    g = RDG.load_example(g)
    assert g.get_edges() == [1,2,3,4]

def test_get_edges_from_to():
    g = RDG()
    g = RDG.load_example(g)
    edges = g.get_edges_from_to()    
    assert sorted(edges) == sorted([(1,3),(3,4), (4,5), (3,2)])