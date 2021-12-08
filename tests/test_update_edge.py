from RDG import RDG, Node, Edge


def test_update_edge_input_to_new_to_node():
    g = RDG()
    g = RDG.load_example(g)
    node_key = g.get_new_node_key()

    start_codon_position = 15
    stop_codon_position = 25
    start_node = Node(key=node_key, node_type="start_codon", coordinates=start_codon_position, edges_in=[], edges_out=[], nodes_in=[], nodes_out=[])

    node_key = node_key + 1
    stop_node = Node(key=node_key, node_type="stop_codon", coordinates=stop_codon_position, edges_in=[], edges_out=[], nodes_in=[], nodes_out=[])
    
    g.update_edge(2, 1, 5, (1,1000))
    tf = 2 in g.edges
    assert 2 in g.nodes[5].input_edges

def test_update_edge_input_to_old_to_node():
    g = RDG()
    g = RDG.load_example(g)
    g.update_edge(2, 1, 5, (1,1000))
    assert 2 not in g.nodes[2].input_edges


def test_update_edge_input_to_new_from_node():
    g = RDG()
    g = RDG.load_example(g)

    g.update_edge(2, 1, 5, (1,1000))
    assert 2 in g.nodes[1].output_edges


def test_update_edge_input_to_old_from_node():
    g = RDG()
    g = RDG.load_example(g)
    g.update_edge(2, 1, 5, (1,1000))
    assert 2 not in g.nodes[3].output_edges


def test_update_edge_input_to_new_to_node():
    g = RDG()
    g = RDG.load_example(g)

    g.update_edge(2, 1, 5, (1,1000))
    assert 2 in g.nodes[1].output_edges