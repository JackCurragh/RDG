"""
This script contains functions for dynamically builidng RDGs 

RDGs can be built 'dynamically' by adding nodes and edges to the graph for the inclusion of features 
or 'statically' where the structure of the graph is pre-defined by the annotation of branch points.  
"""

from RDG import RDG, Node, Edge


def insert_ORF(self, edge: Edge, start_node: Node, stop_node: Node):
    """
    Handle inserting node into just one edge including down path

    Parameters:
    -----------
    edge: Edge object to insert into
    start_node: Node object of the start codon for the ORF
    stop_node: Node object of the stop codon for the ORF

    Notes:
    ------
    Results in the creation of three new edges and two new nodes for the coding branch
    The non-coding branch is updated to now start at the start codon

    """
    five_prime_edge_coords = (edge.coordinates[0], start_node.node_start - 1)
    coding_edge_coords = (start_node.node_start, stop_node.node_start)
    three_prime_edge_coords = (stop_node.node_start + 1, edge.coordinates[1])

    edge_key = self.get_new_edge_key()
    five_prime = Edge(
        key=edge_key,
        edge_type="untranslated",
        from_node=edge.from_node,
        to_node=start_node.key,
        coordinates=five_prime_edge_coords,
    )

    self.add_node(start_node)
    self.add_edge(five_prime, edge.from_node, start_node.key)

    edge_key = self.get_new_edge_key()
    coding = Edge(
        key=edge_key,
        edge_type="translated",
        from_node=start_node.key,
        to_node=stop_node.key,
        coordinates=coding_edge_coords,
    )

    self.add_node(stop_node)
    self.add_edge(coding, start_node.key, stop_node.key)

    terminal_node_key = self.get_new_node_key()
    terminal_node = Node(
        key=terminal_node_key,
        node_type="3_prime",
        position=edge.coordinates[1],
        edges_in=[],
        edges_out=[],
        nodes_in=[],
        nodes_out=[],
    )

    edge_key = self.get_new_edge_key()
    three_prime = Edge(
        key=edge_key,
        edge_type="untranslated",
        from_node=stop_node.key,
        to_node=terminal_node.key,
        coordinates=three_prime_edge_coords,
    )

    self.add_node(terminal_node)
    self.add_edge(three_prime, stop_node.key, terminal_node.key)

    self.update_edge(
        edge.key,
        start_node.key,
        edge.to_node,
        (start_node.node_start, self.nodes[edge.to_node].node_start),
    )


def add_open_reading_frame(
    graph,
    start_codon_position: int,
    stop_codon_position: int,
    reinitiation: bool = False,
    upstream_limit: int = 0,
):
    """
    Handles all operations related to adding a new decision to the graph. Includes objects in graph and corrects all objects involved
    First finds edges that clash with the new ORF (i.e the start codon is within the range of the edge)
    Clashes are only considered if the edge is not of type 'translated'
    For each clash it will insert an ORF. To do so it will create:
        a new start node: the start codon position
        a new stop node: the stop codon position
        a new terminal node: the 3' end of the path
        a new coding edge: from the start node to the stop node
        2 new 3' edges. One from the stop node to the new terminal node and one from the start node to the old terminal node

    Parameters:
    -----------
    start_codon_position: int position of the start codon
    stop_codon_position: int position of the stop codon
    reinitiation: bool whether or not this is allowed to be a reinitiation event
    upstream_limit: int number of ORFs allowed upstream

    """
    exisiting_edges = graph.get_edges()
    clashing_edges = []
    for edge in exisiting_edges:
        if (
            start_codon_position
            in range(graph.edges[edge].coordinates[0], graph.edges[edge].coordinates[1])
            and graph.edges[edge].edge_type != "translated"
        ):
            upstream_node = graph.edges[edge].from_node
            clashing_edges.append((edge, upstream_node))

    for edge, upstream_node in clashing_edges:
        node_key = graph.get_new_node_key()
        start_node = Node(
            key=node_key,
            node_type="start",
            position=start_codon_position,
            edges_in=[],
            edges_out=[],
            nodes_in=[],
            nodes_out=[],
        )

        node_key = node_key + 1
        stop_node = Node(
            key=node_key,
            node_type="stop",
            position=stop_codon_position,
            edges_in=[],
            edges_out=[],
            nodes_in=[],
            nodes_out=[],
        )
        if reinitiation or not graph.check_translation_upstream(
            upstream_node, upstream_limit=upstream_limit
        ):
            graph.insert_ORF(graph.edges[edge], start_node, stop_node)


def add_stop_codon_readthrough(
    graph, readthrough_codon_position: int, next_stop_codon_position: int
):
    """
    Handles all operations related to adding a new decision to the graph at a stop codon. Includes objects in graph and corrects all objects involved

    Here the stop node at the specified position becomes a branch point and a new stop node is added downstream. The edges eminating from the old stop node
    go to new stop node and a old terminal node.

    Parameters:
    -----------
    readthrough_codon_position: int position of the stop codon
    next_stop_codon_position: int position of the next stop codon

    """
    readthrough_codon_keys = graph.get_key_from_position(
        readthrough_codon_position, "stop"
    )

    for readthrough_key in readthrough_codon_keys:
        # Handle coding edge from old stop to new stop
        # old stop is readthrough_key
        # new stop is new_stop_node_key
        new_stop_node_key = graph.get_new_node_key()

        coding_edge_key = graph.get_new_edge_key()

        new_stop_node = Node(
            key=new_stop_node_key,
            node_type="stop",
            position=next_stop_codon_position,
        )
        graph.add_node(new_stop_node)

        coding = Edge(
            key=coding_edge_key,
            edge_type="translated",
            from_node=readthrough_key,
            to_node=new_stop_node.key,
            coordinates=(
                graph.nodes[readthrough_key].node_start,
                next_stop_codon_position,
            ),
        )
        graph.add_edge(coding, readthrough_key, new_stop_node.key)

        graph.nodes[readthrough_key].node_type = "stop"

        terminal_node_key = graph.get_new_node_key()
        three_prime_terminal_key = graph.nodes[readthrough_key].output_nodes[0]

        new_3_prime_edge = graph.get_new_edge_key()
        three_prime = Edge(
            key=new_3_prime_edge,
            edge_type="untranslated",
            from_node=new_stop_node.key,
            to_node=terminal_node_key,
            coordinates=(
                new_stop_node.node_start,
                graph.nodes[three_prime_terminal_key].node_start,
            ),
        )

        terminal_node_key = graph.get_new_node_key()
        terminal_node = Node(
            key=terminal_node_key,
            node_type="3_prime",
            position=graph.nodes[three_prime_terminal_key].node_start,
        )
        graph.add_node(terminal_node)
        graph.add_edge(three_prime, new_stop_node.key, terminal_node_key)

        # terminal_node_key2 = graph.get_new_node_key()

        # terminal_node2 = Node(
        #     key=terminal_node_key2,
        #     node_type="3_prime",
        #     position=graph.nodes[three_prime_terminal_key].node_start,
        #     edges_in=[new_3_prime_edge],
        #     edges_out=[],
        #     nodes_in=[new_stop_node.key],
        #     nodes_out=[],
        # )
        # graph.add_node(terminal_node2)
        # graph.add_edge(three_prime, readthrough_key, terminal_node_key2)

        # print(graph.nodes)


def add_frameshift(graph, fs_position, next_stop_codon_position, shift):
    """
    add a frameshifting event
    """
    exisiting_edges = graph.get_edges()
    clashing_edges = []
    for edge in exisiting_edges:
        if (
            fs_position
            in range(graph.edges[edge].coordinates[0], graph.edges[edge].coordinates[1])
            and graph.edges[edge].edge_type == "translated"
        ):
            upstream_node = graph.edges[edge].from_node
            clashing_edges.append((edge, upstream_node))

    for edge, upstream_node in clashing_edges:
        # edge and upstrream are keys
        shift_node_key = graph.get_new_node_key()
        new_stop_node_key = shift_node_key + 1
        old_stop_node_key = graph.edges[edge].to_node
        new_terminal_key = new_stop_node_key + 1

        old_stop_edge_key = graph.get_new_edge_key()
        new_stop_edge_key = old_stop_edge_key + 1
        new_three_prime_edge_key = new_stop_edge_key + 1

        shift_node = Node(
            key=shift_node_key,
            node_type="frameshift",
            coordinates=fs_position,
            edges_in=[edge],
            edges_out=[old_stop_edge_key, new_stop_edge_key],
            nodes_in=[upstream_node],
            nodes_out=[old_stop_node_key, new_stop_node_key],
        )
        graph.add_node(shift_node)

        graph.nodes[upstream_node].output_nodes.remove(old_stop_node_key)
        graph.nodes[upstream_node].output_nodes.append(shift_node.key)

        graph.edges[edge].to_node = shift_node_key

        graph.nodes[old_stop_node_key].input_edges.remove(edge)
        graph.nodes[old_stop_node_key].input_edges.append(old_stop_edge_key)
        graph.nodes[old_stop_node_key].input_nodes.remove(upstream_node)
        graph.nodes[old_stop_node_key].input_nodes.append(shift_node_key)

        new_stop_node = Node(
            key=new_stop_node_key,
            node_type="stop",
            coordinates=next_stop_codon_position,
            edges_in=[new_stop_edge_key],
            edges_out=[new_three_prime_edge_key],
            nodes_in=[shift_node_key],
            nodes_out=[new_terminal_key],
        )
        graph.add_node(new_stop_node)

        new_terminal_node = Node(
            key=new_terminal_key,
            node_type="3_prime",
            coordinates=graph.locus_stop,
            edges_in=[new_three_prime_edge_key],
            edges_out=[],
            nodes_in=[new_stop_node_key],
            nodes_out=[],
        )
        graph.add_node(new_terminal_node)

        old_stop_edge = Edge(
            key=old_stop_edge_key,
            edge_type="translated",
            from_node=shift_node_key,
            to_node=old_stop_node_key,
            coordinates=(
                shift_node.node_start,
                graph.nodes[old_stop_node_key].node_start,
            ),
        )
        graph.add_edge(
            old_stop_edge,
            from_node_key=shift_node_key,
            to_node_key=old_stop_node_key,
        )
        graph.edges[old_stop_edge_key].frame = graph.edges[edge].frame

        new_stop_edge = Edge(
            key=new_stop_edge_key,
            edge_type="translated",
            from_node=shift_node_key,
            to_node=new_stop_node_key,
            coordinates=(
                shift_node.node_start,
                graph.nodes[new_stop_node_key].node_start,
            ),
        )
        graph.add_edge(
            new_stop_edge,
            from_node_key=shift_node_key,
            to_node_key=new_stop_node_key,
        )

        new_three_prime_edge = Edge(
            key=new_three_prime_edge_key,
            edge_type="untranslated",
            from_node=new_stop_node_key,
            to_node=new_terminal_key,
            coordinates=(new_stop_node.node_start, new_terminal_node.node_start),
        )
        graph.add_edge(
            new_three_prime_edge,
            from_node_key=new_stop_node_key,
            to_node_key=new_terminal_key,
        )
