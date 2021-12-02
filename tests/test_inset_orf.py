from RDG import RDG, Node, Edge

def test_inserting_orfs_edges():
    g = RDG()
    g = RDG.load_example(g)
    node_key = g.get_new_node_key()

    start_codon_position = 15
    stop_codon_position = 25
    start_node = Node(key=node_key, node_type="start_codon", coordinates=(start_codon_position, start_codon_position + 2), edges_in=[], edges_out=[], nodes_in=[], nodes_out=[])

    node_key = node_key + 1
    stop_node = Node(key=node_key, node_type="stop_codon", coordinates=(stop_codon_position, stop_codon_position + 2), edges_in=[], edges_out=[], nodes_in=[], nodes_out=[])
    
    g.insert_ORF(g.edges[2], start_node, stop_node)
    assert list(g.edges.keys()) == [1,2,3,4,5,6,7]


def test_insterting_orfs_nodes():
    g = RDG()
    g = RDG.load_example(g)
    node_key = g.get_new_node_key()

    start_codon_position = 15
    stop_codon_position = 25
    start_node = Node(key=node_key, node_type="start_codon", coordinates=(start_codon_position, start_codon_position + 2), edges_in=[], edges_out=[], nodes_in=[], nodes_out=[])

    node_key = node_key + 1
    stop_node = Node(key=node_key, node_type="stop_codon", coordinates=(stop_codon_position, stop_codon_position + 2), edges_in=[], edges_out=[], nodes_in=[], nodes_out=[])
    
    g.insert_ORF(g.edges[2], start_node, stop_node)
    assert list(g.nodes.keys()) == [1,2,3,4,5,6,7,8]


def test_insterting_orfs_branch_points():
    g = RDG()
    g = RDG.load_example(g)
    node_key = g.get_new_node_key()

    start_codon_position = 15
    stop_codon_position = 25
    start_node = Node(key=node_key, node_type="start_codon", coordinates=(start_codon_position, start_codon_position + 2), edges_in=[], edges_out=[], nodes_in=[], nodes_out=[])

    node_key = node_key + 1
    stop_node = Node(key=node_key, node_type="stop_codon", coordinates=(stop_codon_position, stop_codon_position + 2), edges_in=[], edges_out=[], nodes_in=[], nodes_out=[])
    
    g.insert_ORF(g.edges[2], start_node, stop_node)
    
    branch_points = g.get_branch_points()

    assert len(branch_points) == 2