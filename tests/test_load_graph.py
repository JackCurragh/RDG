from RDG import RDG, Node, Edge


def test_load_locus_stop_update():
    g = RDG()
    g = RDG.load(g, locus_name="GeneA", locus_start=0, locus_stop=10000)
    assert g.locus_stop == 10000


def test_load_nodes():
    g = RDG()
    g = RDG.load_example(g)
    assert len(list(g.nodes.keys())) == 5