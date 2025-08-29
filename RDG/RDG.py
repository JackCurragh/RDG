"""
This script contains the RDG class and associated classes for nodes and edges.

The RDG class is the main class for the RDG package. It contains methods for
creating, loading, and saving RDGs. It also contains methods for adding and
removing nodes and edges from the graph.

The Node class is a class for nodes in the RDG. It contains methods for
getting information about the node.

The Edge class is a class for edges in the RDG. It contains methods for
getting information about the edge.

This is the first test implementation of the RDG concept and is primarily
intended for testing and development.
"""

from typing import Tuple, Dict, List
from dataclasses import dataclass, field

@dataclass
class Node:
    """
    Represents a node in a decision graph.

    Attributes:
        key (str): A unique identifier for the node.
        node_type (str): Type of the node, specifying its function or role.
        input_edges (list): List of edges incoming to the node.
        output_edges (list): List of edges outgoing from the node.
        input_nodes (list): List of nodes connected to the input edges.
        output_nodes (list): List of nodes connected to the output edges.
        node_start (int): Position of the node within the sequence.
        frame (int): Reading frame of the node within the sequence.
    """
    key: int
    node_type: str
    position: int
    input_edges: List[int] = field(default_factory=list)
    output_edges: List[int] = field(default_factory=list)
    input_nodes: List[int] = field(default_factory=list)
    output_nodes: List[int] = field(default_factory=list)

    def __post_init__(self):
        """Validate the node_type after initialization."""
        valid_types = [
            "5_prime",
            "3_prime",
            "start",
            "stop",
            "frameshift",
            "readthrough_stop",
        ]
        if self.node_type not in valid_types:
            raise ValueError(f"Invalid node type: {self.node_type}")

    @property
    def frame(self) -> int:
        """
        Get the reading frame of the node within the sequence.

        Returns:
            int: The reading frame of the node.
        """
        return self.position % 3


@dataclass
class Edge:
    """
    Represents an edge in a decision graph.

    Attributes:
        key (int): A unique identifier for the edge.
        edge_type (str): Type of the edge, specifying its function or role.
        from_node (int): The source node of the edge.
        to_node (int): The destination node of the edge.
        coordinates (Tuple[int, int]): Tuple representing the start and
                                        stop positions of the edge.
    """
    key: int
    edge_type: str
    from_node: int
    to_node: int
    coordinates: Tuple[int, int]

    @property
    def frame(self) -> int:
        """
        Get the reading frame of the edge within the sequence.

        Returns:
            int: The reading frame of the edge.
        """
        return self.coordinates[0] % 3
    
@dataclass
class RDG:
    """
    Represents a Reading Decision Graph (RDG).
    
    Attributes:
        locus (str): Name of the locus.
        locus_start (int): Start of the locus in base pairs.
        locus_stop (int): Length of the locus in base pairs.
        nodes (Dict[int, Node]): Dictionary of nodes. Keys are node keys, values are Node objects.
        edges (Dict[int, Edge]): Dictionary of edges. Keys are edge keys, values are Edge objects.
    """
    name: str = "name"
    locus_start: int = 0
    locus_stop: int = 1000
    nodes: Dict[int, Node] = field(default_factory=dict)
    edges: Dict[int, Edge] = field(default_factory=dict)

    def __post_init__(self):
        """Initialize a fresh RDG object if no nodes are provided."""
        if self.locus_stop <= self.locus_start:
            raise ValueError("locus_stop must be greater than locus_start")
            
        if not self.nodes:
            # Initialize a new RDG object with default nodes and edges
            self.nodes = {
                1: Node(
                    key=1,
                    node_type="5_prime",
                    position=0,
                    output_edges=[1],
                    output_nodes=[2]
                ),
                2: Node(
                    key=2,
                    node_type="3_prime",
                    position=self.locus_stop - self.locus_start,
                    input_edges=[1],
                    input_nodes=[1],
                ),
            }
            self.edges = {
                1: Edge(
                    key=1,
                    edge_type="untranslated",
                    from_node=1,
                    to_node=2,
                    coordinates=(1, self.locus_stop - 1),
                )
            }
        self._position_index = {}
        self._rebuild_position_index()

    def _rebuild_position_index(self):
        """Rebuild the position index for the RDG."""
        self._position_index = {}
        for node_key, node in self.nodes.items():
            key = (node.position, node.node_type)
            if key not in self._position_index:
                self._position_index[key] = []
            self._position_index[key].append(node_key)


    def load_example(self) -> 'RDG':
        """
        Load a basic graph with one translon.

        Returns:
        RDG object: RDG object with the following attributes:
            - locus length = 1000
            - translon coordinates = (10, 100)
        """
        self.nodes = {
            1: Node(
                key=1,
                node_type="5_prime",
                position=0,
                output_edges=[1],
                output_nodes=[3]
                ),
            2: Node(
                key=2,
                node_type="3_prime",
                position=1000,
                input_edges=[2],
                input_nodes=[3]
                ),
            3: Node(
                key=3,
                node_type="start",
                position=10,
                input_edges=[1],
                output_edges=[2, 3],
                input_nodes=[1],
                output_nodes=[2, 4],
            ),
            4: Node(
                key=4,
                node_type="stop",
                position=100,
                input_edges=[3],
                output_edges=[4],
                input_nodes=[3],
                output_nodes=[5],
            ),
            5: Node(
                key=5,
                node_type="3_prime",
                position=1000,
                input_edges=[4],
                input_nodes=[4]
                ),
        }

        self.edges = {
            1: Edge(
                1,
                "untranslated",
                from_node=1,
                to_node=3,
                coordinates=(1, 10 - 1)
                ),
            2: Edge(
                2,
                "untranslated",
                from_node=3,
                to_node=2,
                coordinates=(11, 1000 - 1)
                ),
            3: Edge(
                3,
                "translated",
                from_node=3,
                to_node=4,
                coordinates=(11, 100)
                ),
            4: Edge(
                4,
                "untranslated",
                from_node=4,
                to_node=5,
                coordinates=(101, 1000 - 1)
                ),
        }

        self._rebuild_position_index()
        return self

    def get_node_keys(self) -> List[int]:
        """
        Get the keys of all nodes in the graph.

        Returns:
        List[int]: A list containing the keys of all nodes.
        """
        return list(self.nodes.keys())

    def get_edge_keys(self) -> List[int]:
        """
        Get the keys of all edges in the graph.

        Returns:
        List[int]: A list containing the keys of all edges.
        """
        return list(self.edges.keys())

    def get_edges_from_to(self) -> dict:
        """
        Return a dict of edges with keys of the form (from_node_id, to_node_id)

        Returns:
        dict Keys: tuples: (from_node_id, to_node_id) Values: edge ids
        """
        return {
            (self.edges[edge].from_node, self.edges[edge].to_node): edge
            for edge in self.edges.keys()
        }

    def get_new_node_key(self) -> int:
        """
        Return the lowest available node key

        Returns:
        int: new node key that is not already used
        """
        keys = self.nodes.keys()
        if len(keys) > 0:
            return max(keys) + 1
        else:
            return 1

    def get_new_edge_key(self) -> int:
        """
        Return the lowest available edge key

        Returns:
        int: new edge key that is not already used
        """
        keys = self.edges.keys()
        if len(keys) > 0:
            return max(keys) + 1
        else:
            return 1

    def get_key_from_position(self, position: int, node_type: str) -> list:
        """
        Return the node key for the node of specified type at the
        stated position.

        Parameters:
        position (int): Position in the locus.
        node_type (str): Type of node to look for.

        Returns:
        int: Node key.
        """
        # Check if any nodes of this type exist at all
        has_type = any((pos, ntype) for pos, ntype in self._position_index.keys() if ntype == node_type)
        if not has_type:
            raise ValueError(f"There are no nodes of type '{node_type}' in the graph")

        key = (position, node_type)
        nodes = self._position_index.get(key, [])
        if not nodes:
            raise ValueError(f"No node of type '{node_type}' at position {position}")
        return nodes

    def remove_edge(self, edge_key: int):
        """
        Remove the edge with the given key from the graph, also removing
        references from adjacent nodes.

        Does not remove nodes that are no longer connected to any other nodes.

        Parameters:
        edge_key (int): Key of the edge to remove.
        """
        from_node = self.edges[edge_key].from_node
        to_node = self.edges[edge_key].to_node
        if edge_key in self.nodes[from_node].output_edges:
            self.nodes[from_node].output_edges.remove(edge_key)

        if to_node in self.nodes[from_node].output_nodes:
            self.nodes[from_node].output_nodes.remove(to_node)

        if edge_key in self.nodes[to_node].input_edges:
            self.nodes[to_node].input_edges.remove(edge_key)

        if from_node in self.nodes[to_node].input_nodes:
            self.nodes[to_node].input_nodes.remove(from_node)

        if edge_key in self.edges:
            self.edges.pop(edge_key)

    def add_edge(self, edge: Edge, from_node_key: int, to_node_key: int):
        """
        Introduce an edge and its references to the graph.

        Parameters:
        edge (Edge): Edge object to add.
        from_node_key (int): Key of the node the edge is coming from.
        to_node_key (int): Key of the node the edge is going to.
        """

        self.nodes[from_node_key].output_edges.append(edge.key)
        self.nodes[from_node_key].output_nodes.append(to_node_key)

        self.nodes[to_node_key].input_edges.append(edge.key)
        self.nodes[to_node_key].input_nodes.append(from_node_key)

        self.edges[edge.key] = edge

    def add_node(self, node: Node):
        """
        Introduce nodes and references with the a node.

        The associated edges must be added separately

        Parameters:
        node: Node object to add

        """
        self.nodes[node.key] = node

        key = (node.position, node.node_type)
        if key not in self._position_index:
            self._position_index[key] = []
        self._position_index[key].append(node.key)

    def remove_node(self, node_key: int):
        """
        Remove the node with the given key from the graph also removing
        references from adjacent edges

        Does not remove edges that are no longer connected to any other nodes

        Parameters:
        node_key: int key of the node to remove
        """
        if node_key in self.nodes:
            node = self.nodes[node_key]
            key = (node.position, node.node_type)
            if key in self._position_index:
                self._position_index[key].remove(node_key)
                if not self._position_index[key]:
                    del self._position_index[key]

            for edge in sorted(self.nodes[node_key].output_edges):
                self.remove_edge(edge)
            for edge in sorted(self.nodes[node_key].input_edges):
                self.remove_edge(edge)
            self.nodes.pop(node_key)

    def update_edge(
            self,
            edge_key: int,
            new_from_node: int,
            new_to_node: int,
            new_coordinates: Tuple[int, int]
            ):
        """
        Update the edge with the given key with the new from and to
        nodes and coordinates
        notes:
        when inserting an translon into an untranslated edge we update the
        untranslated edge to now start at the new start codon instead

        Parameters:
        edge_key: int key of the edge to update
        new_from_node: int key of the new from node
        new_to_node: int key of the new to node
        new_coordinates: tuple of the form (start, end) of the new coordinates
        """
        old_from_node = self.edges[edge_key].from_node
        old_to_node = self.edges[edge_key].to_node

        if old_from_node != new_from_node:
            if edge_key in self.nodes[old_from_node].output_edges:
                self.nodes[old_from_node].output_edges.remove(edge_key)

            if edge_key not in self.nodes[new_from_node].output_edges:
                self.nodes[new_from_node].output_edges.append(edge_key)

            if old_to_node in self.nodes[old_from_node].output_nodes:
                self.nodes[old_from_node].output_nodes.remove(old_to_node)

            if new_from_node not in self.nodes[new_from_node].output_nodes:
                self.nodes[new_from_node].output_nodes.append(new_to_node)

            if old_from_node in self.nodes[new_to_node].input_nodes:
                self.nodes[new_to_node].input_nodes.remove(old_from_node)

            if new_from_node not in self.nodes[new_to_node].input_nodes:
                self.nodes[new_to_node].input_nodes.append(new_from_node)

        if old_to_node != new_to_node:
            if edge_key in self.nodes[old_to_node].input_edges:
                self.nodes[old_to_node].input_edges.remove(edge_key)

            if edge_key not in self.nodes[new_to_node].input_edges:
                self.nodes[new_to_node].input_edges.append(edge_key)

            if old_from_node in self.nodes[old_to_node].input_nodes:
                self.nodes[old_to_node].input_nodes.remove(old_from_node)

        self.edges[edge_key].coordinates = new_coordinates
        self.edges[edge_key].from_node = new_from_node

    def insert_translon(self, edge: Edge, start_node: Node, stop_node: Node):
        """
        Handle inserting node into just one edge including down path.
        Nodes are already inserted in the graph but are not connected. This
        function connects them and updates the edges and nodes as necessary
        Results in the creation of three new edges and two new nodes for the
        coding branch. The non-coding branch is updated to now start at the
        start codon

        Parameters:
        edge: Edge object to insert into
        start_node: Node object of the start codon for the translon
        stop_node: Node object of the stop codon for the translon
        """
        five_prime_edge_coords = (
            edge.coordinates[0], start_node.position - 1
            )
        coding_edge_coords = (
            start_node.position, stop_node.position
            )
        three_prime_edge_coords = (
            stop_node.position + 1, self.locus_stop
            )

        edge_key = self.get_new_edge_key()
        five_prime = Edge(
            key=edge_key,
            edge_type="untranslated",
            from_node=edge.from_node,
            to_node=start_node.key,
            coordinates=five_prime_edge_coords,
        )

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
            position=self.locus_stop,
            input_edges=[],
            output_edges=[],
            input_nodes=[],
            output_nodes=[],
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
            (start_node.position, self.nodes[edge.to_node].position),
        )

    def is_input_edge_translated(self, node: int) -> bool:
        """
        Return boolean as to whether an edge entering this node is translated

        Parameters:
        node: int key of the node to check

        Returns:
        boolean: True if an edge entering this node is translated,
                False otherwise
        """
        input_edges = self.nodes[node].input_edges
        boolean = False
        for edge in input_edges:
            if self.edges[edge].edge_type == "translated":
                boolean = True
        return boolean

    def root_to_node_of_acyclic_node_path(self, node: int) -> list:
        """
        Return path from root to given node as a list of node keys

        Parameters:
        node: int key of the node to check

        Returns:
        path: list of node keys from root to given node
        """
        path = [node]
        reached_root = False

        while not reached_root:
            in_node = self.nodes[node].input_nodes[0]
            path.append(in_node)
            if self.nodes[in_node].input_nodes == []:
                reached_root = True
            else:
                node = in_node
        return path[::-1]

    def root_to_node_of_acyclic_edge_path(self, node: int):
        """
        Return path from root to given node as a list of edge keys

        Parameters:
        node: int key of the node to check

        Returns:
        path: list of edge keys from root to given node

        """
        if self.nodes[node].input_nodes == []:
            reached_root = True
            return []

        path = [self.nodes[node].input_edges[0]]
        reached_root = False
        while not reached_root:
            in_node = self.nodes[node].input_nodes[0]
            if self.nodes[in_node].input_nodes == []:
                reached_root = True
            else:
                path.append(self.nodes[in_node].input_edges[0])
                node = in_node
        return path[::-1]

    def check_translation_upstream(
        self, from_node: int, upstream_limit: int = 1
    ) -> bool:
        """
        look upstream of an edges from node and see if any edges are of type
        'translated' used in reinitiation functionality. Upstream limit refers
        to the number of translons allowed upstream. ie. can the be multiple
        reinitiation events or just one etc

        Parameters:
        -----------
        from_node: int key of the node to check
        upstream_limit: int number of translons allowed upstream

        Returns:
        --------
        boolean: True if the number of translated edges upstream is greater
                than the upstream limit, False otherwise
        """
        edge_path_to_node = self.root_to_node_of_acyclic_edge_path(from_node)
        number_of_translated_regions = 0
        for edge in edge_path_to_node:
            if self.edges[edge].edge_type == "translated":
                number_of_translated_regions += 1

        if number_of_translated_regions <= upstream_limit:
            return False
        return True

    def add_open_reading_frame(
        self,
        start_codon_position: int,
        stop_codon_position: int,
        reinitiation: bool = False,
        upstream_limit: int = 0,
    ):
        """
        Handles all operations related to adding a new decision to the graph.
        Includes objects in graph and corrects all objects involved
        First finds edges that clash with the new translon
        (i.e the start codon is within the range of the edge)
        Clashes are only considered if the edge is not of type 'translated'
        For each clash it will insert an translon. To do so it will create:
            a new start node: the start codon position
            a new stop node: the stop codon position
            a new terminal node: the 3' end of the path
            a new coding edge: from the start node to the stop node
            2 new 3' edges. One from the stop node to the new terminal node
            and one from the start node to the old terminal node

        Parameters:
        start_codon_position: int position of the start codon
        stop_codon_position: int position of the stop codon
        reinitiation: bool whether or not this is allowed to be a
            reinitiation event
        upstream_limit: int number of translons allowed upstream

        """
        if stop_codon_position > self.locus_stop:
            raise Exception(
                f"""
                Next in frame stop codon ({stop_codon_position}) is outside
                of sequence (length {self.locus_stop})
                """
            )

        existing_edges = self.get_edge_keys()
        clashing_edges = []
        for edge in existing_edges:
            edge_obj = self.edges[edge]
            if (
                start_codon_position
                in range(
                    edge_obj.coordinates[0], edge_obj.coordinates[1]
                )
                and self.edges[edge].edge_type != "translated"
            ):
                upstream_node = self.edges[edge].from_node
                clashing_edges.append((edge, upstream_node))

        for edge, upstream_node in clashing_edges:
            if reinitiation or not self.check_translation_upstream(
                upstream_node, upstream_limit=upstream_limit
            ):
                node_key = self.get_new_node_key()
                start_node = Node(
                    key=node_key,
                    node_type="start",
                    position=start_codon_position,
                    input_edges=[],
                    output_edges=[],
                    input_nodes=[],
                    output_nodes=[],
                )
                self.add_node(start_node)

                stop_node_key = self.get_new_node_key()
                stop_node = Node(
                    key=stop_node_key,
                    node_type="stop",
                    position=stop_codon_position,
                    input_edges=[],
                    output_edges=[],
                    input_nodes=[],
                    output_nodes=[],
                )
                self.add_node(stop_node)

                self.insert_translon(self.edges[edge], start_node, stop_node)

    def add_stop_codon_readthrough(
        self, readthrough_codon_position: int, next_stop_codon_position: int
    ):
        """
        Handles all operations related to adding a new decision to the graph
        at a stop codon. Includes objects in graph and corrects all objects
        involved

        Here the stop node at the specified position becomes a branch point
        and a new stop node is added downstream. The edges emanating from the
        old stop node go to new stop node and a old terminal node.

        Parameters:
        readthrough_codon_position: int position of the stop codon
        next_stop_codon_position: int position of the next stop codon

        """
        if next_stop_codon_position > self.locus_stop:
            raise Exception(
                f"""
                Next in frame stop codon ({next_stop_codon_position}) is
                outside of sequence (length {self.locus_stop})
                """
            )

        readthrough_codon_keys = self.get_key_from_position(
            readthrough_codon_position, "stop"
        )

        for readthrough_key in readthrough_codon_keys:
            # Handle coding edge from old stop to new stop
            # old stop is readthrough_key
            # new stop is new_stop_node_key
            new_stop_node_key = self.get_new_node_key()
            new_stop_node = Node(
                key=new_stop_node_key,
                node_type="stop",
                position=next_stop_codon_position,
                input_edges=[],
                output_edges=[],
                input_nodes=[],
                output_nodes=[],
            )
            self.add_node(new_stop_node)

            coding_edge_key = self.get_new_edge_key()
            coding = Edge(
                key=coding_edge_key,
                edge_type="translated",
                from_node=readthrough_key,
                to_node=new_stop_node.key,
                coordinates=(
                    self.nodes[readthrough_key].position,
                    next_stop_codon_position,
                ),
            )
            self.add_edge(coding, readthrough_key, new_stop_node.key)

            # Handle 3' edge from new stop to new terminal node
            terminal_node_key = self.get_new_node_key()
            three_prime_terminal_key = self.nodes[
                readthrough_key
                ].output_nodes[0]

            terminal_node_key = self.get_new_node_key()
            terminal_node = Node(
                key=terminal_node_key,
                node_type="3_prime",
                position=self.nodes[three_prime_terminal_key].position,
                input_nodes=[],
                output_nodes=[],
                input_edges=[],
                output_edges=[],
            )
            self.add_node(terminal_node)

            new_3_prime_edge = self.get_new_edge_key()
            three_prime = Edge(
                key=new_3_prime_edge,
                edge_type="untranslated",
                from_node=new_stop_node.key,
                to_node=terminal_node_key,
                coordinates=(
                    new_stop_node.position,
                    self.nodes[three_prime_terminal_key].position,
                ),
            )
            self.add_edge(three_prime, new_stop_node.key, terminal_node_key)

    def add_frameshift(
            self,
            fs_position: int,
            next_stop_codon_position: int,
            shift: int
            ) -> None:
        """
        Handles adding a frameshifting event to the graph.

        Frameshifts can only be added to translated edges. If the frameshift
        is added to a translated edge, the edge is split into two edges. The
        first edge is the edge up from the frameshift position to the stop.
        The second edge is the original edge from the frameshift position
        to the stop codon. The first edge is then translated and the second
        edge is untranslated.

        Parameters:
        fs_position: int position of the frameshift
        next_stop_codon_position: int position of the next stop codon after
        the frameshift
        shift: int number of bases to shift by (-1, +1 etc.)
        """

        if next_stop_codon_position > self.locus_stop:
            raise Exception(
                f"""
                Next in frame stop codon ({next_stop_codon_position}) is
                outside of sequence (length {self.locus_stop})
                """
            )

        existing_edges = self.get_edge_keys()
        clashing_edges = []
        for edge in existing_edges:
            edge_obj = self.edges[edge]
            if (
                fs_position
                in range(
                    edge_obj.coordinates[0], edge_obj.coordinates[1]
                )
                and edge_obj.edge_type == "translated"
            ):
                upstream_node = edge_obj.from_node
                clashing_edges.append((edge, upstream_node))

        for edge, upstream_node in clashing_edges:
            # Add FS node
            shift_node_key = self.get_new_node_key()
            shift_node = Node(
                key=shift_node_key,
                node_type="frameshift",
                position=fs_position,
                input_edges=[],
                output_edges=[],
                input_nodes=[],
                output_nodes=[],
            )
            self.add_node(shift_node)

            # Add FS edge to old stop node (i.e the event where no FS happened)
            old_stop_node_key = self.edges[edge].to_node

            old_stop_edge_key = self.get_new_edge_key()
            old_stop_edge = Edge(
                key=old_stop_edge_key,
                edge_type="translated",
                from_node=shift_node_key,
                to_node=old_stop_node_key,
                coordinates=(
                    shift_node.position + shift,
                    self.nodes[old_stop_node_key].position,
                ),
            )
            self.add_edge(
                old_stop_edge,
                from_node_key=shift_node_key,
                to_node_key=old_stop_node_key,
            )

            # update the upstream node and old stop node so they have correct
            # references reflecting the addition of a FS
            self.nodes[upstream_node].output_nodes.remove(old_stop_node_key)
            self.nodes[upstream_node].output_nodes.append(shift_node.key)

            self.edges[edge].to_node = shift_node_key

            self.nodes[old_stop_node_key].input_edges.remove(edge)
            self.nodes[old_stop_node_key].input_edges.append(old_stop_edge_key)

            if upstream_node in self.nodes[old_stop_node_key].input_nodes:
                self.nodes[old_stop_node_key].input_nodes.remove(upstream_node)

            # Add FS edge to new stop node (i.e the event where a FS happened)
            new_stop_node_key = self.get_new_node_key()
            new_stop_node = Node(
                key=new_stop_node_key,
                node_type="stop",
                position=next_stop_codon_position,
                input_edges=[],
                output_edges=[],
                input_nodes=[],
                output_nodes=[],
            )
            self.add_node(new_stop_node)

            new_stop_edge_key = self.get_new_edge_key()
            new_stop_edge = Edge(
                key=new_stop_edge_key,
                edge_type="translated",
                from_node=shift_node_key,
                to_node=new_stop_node_key,
                coordinates=(
                    shift_node.position,
                    self.nodes[new_stop_node_key].position,
                ),
            )
            self.add_edge(
                new_stop_edge,
                from_node_key=shift_node_key,
                to_node_key=new_stop_node_key,
            )

            # Add new terminal node and edge to that node
            # from the new stop node
            new_terminal_key = self.get_new_node_key()
            new_terminal_node = Node(
                key=new_terminal_key,
                node_type="3_prime",
                position=self.locus_stop,
                input_edges=[],
                output_edges=[],
                input_nodes=[],
                output_nodes=[],
            )
            self.add_node(new_terminal_node)

            new_three_prime_edge_key = self.get_new_edge_key()
            new_three_prime_edge = Edge(
                key=new_three_prime_edge_key,
                edge_type="untranslated",
                from_node=new_stop_node_key,
                to_node=new_terminal_key,
                coordinates=(
                    new_stop_node.position,
                    new_terminal_node.position
                    ),
            )
            self.add_edge(
                new_three_prime_edge,
                from_node_key=new_stop_node_key,
                to_node_key=new_terminal_key,
            )

    def get_branch_points(self) -> list:
        """
        Return a list of nodes that are branch points
        (i.e have more than 1 output edge)

        Returns:
        list
        """
        branch_points = []
        for node in self.nodes:
            if len(self.nodes[node].output_edges) > 1:
                branch_points.append(node)
        return branch_points

    def get_endpoints(self) -> list:
        """
        Return list of terminal node keys. These nodes have no output edges as
        they are the end of the path through the graph

        Returns:
        list
        """
        endpoints = []
        for node in self.nodes:
            if (
                len(self.nodes[node].output_edges) == 0
                and self.nodes[node].node_type == "3_prime"
            ):
                endpoints.append(node)
        return endpoints

    def get_startpoints(self) -> list:
        """
        Return list of start node keys. These nodes have no input edges as
        they are the start of the path through the graph (i.e the 5' end of
        the locus or TSS)

        Returns:
        list
        """
        startpoints = []
        for node in self.nodes:
            if (
                len(self.nodes[node].input_edges) == 0
                and self.nodes[node].node_type == "5_prime"
            ):
                startpoints.append(node)
        return startpoints

    def get_start_nodes(self) -> list:
        """
        Return translation start keys. These nodes have 2 output edges. One
        for where translation starts at that site and one for where it doesn't

        Returns:
        list
        """
        translation_starts = []
        for node in self.nodes:
            if (
                len(self.nodes[node].output_edges) == 2
                and self.nodes[node].node_type == "start"
            ):
                translation_starts.append(node)
        return translation_starts

    def get_stop_nodes(self) -> list:
        """
        Return translation stop keys. These may have 1 or 2 output edges.
        If they have 2 output edges then they are readthrough cases

        Returns:
        list
        """
        translation_stops = []
        for node in self.nodes:
            if "stop" in self.nodes[node].node_type:
                translation_stops.append(node)
        return translation_stops

    def get_translons(self) -> list:
        """
        Return a list of coordinates (start, stop) describing the translons in
        the graph. Frameshift translons have the form (start, FS coordinate,
        stop) even when FS event is negative

        NOTE: multiple readthroughs are not yet supported

        Returns:
        list
        """
        translation_starts = self.get_start_nodes()

        translons = []
        for start in translation_starts:
            downstream_nodes = self.nodes[start].output_nodes

            for candidate_stop in downstream_nodes:
                if self.nodes[candidate_stop].node_type == "stop":
                    translons.append(
                        (
                            self.nodes[start].position,
                            self.nodes[candidate_stop].position,
                        )
                    )

                    # Search for cases of readthrough
                    readthrough_downstream_nodes = self.nodes[
                        candidate_stop
                    ].output_nodes

                    for new_candidate_stop in readthrough_downstream_nodes:
                        if self.nodes[new_candidate_stop].node_type == "stop":
                            translon = (
                                self.nodes[start].position,
                                self.nodes[new_candidate_stop].position,
                            )
                            if translon not in translons:
                                translons.append(translon)

                elif self.nodes[candidate_stop].node_type == "frameshift":
                    fs_downstream_nodes = self.nodes[
                        candidate_stop
                        ].output_nodes

                    for new_candidate_stop in fs_downstream_nodes:
                        if self.nodes[new_candidate_stop].node_type == "stop":
                            translon = (
                                self.nodes[start].position,
                                self.nodes[candidate_stop].position,
                                self.nodes[new_candidate_stop].position,
                            )
                            if translon not in translons:
                                translons.append(translon)

        return translons

    def get_frameshifts(self) -> list:
        """
        Return a list of keys for all frameshift nodes.

        Returns:
        list
        """
        frameshifts = []
        for node in self.nodes:
            if "frameshift" in self.nodes[node].node_type:
                frameshifts.append(node)
        return frameshifts

    def get_unique_paths(self) -> list:
        """
        Return a list of lists of nodes that describe the unique
        paths through the graph

        Returns:
        list
        """
        terminal_nodes = self.get_endpoints()
        paths = []
        for node in terminal_nodes:
            path_to_root = self.root_to_node_of_acyclic_node_path(node)
            if path_to_root not in paths:
                paths.append(path_to_root)
        return paths

    def statistics(self) -> dict:
        """
        Return a dictionary of statistics about the graph

        Returns:
        dict
        """
        stats = {}
        nodes = list(self.nodes.keys())
        calc_dict = {
            "edges": [],
            "frames": [],
            "types": [],
            "frames_freq": {},
            "types_freq": {},
        }

        calc_dict["edges"] = self.get_edges_from_to().keys()

        for i in nodes:
            calc_dict["frames"].append(self.nodes[i].frame)
            calc_dict["types"].append(self.nodes[i].node_type)

        for item in calc_dict["frames"]:
            if item in calc_dict["frames_freq"]:
                calc_dict["frames_freq"][item] += 1
            else:
                calc_dict["frames_freq"][item] = 1

        for item in calc_dict["types"]:
            if item in calc_dict["types_freq"]:
                calc_dict["types_freq"][item] += 1
            else:
                calc_dict["types_freq"][item] = 1

        for item in calc_dict["types_freq"]:
            key = "Number_of_" + item
            stats[key] = calc_dict["types_freq"][item]

        for item in calc_dict["frames_freq"]:
            key = "Number_of_input_nodes_frame_" + str(item)
            stats[key] = calc_dict["frames_freq"][item]

        stats["Node_keys"] = nodes
        stats["Number_of_nodes"] = len(self.nodes.keys())
        stats["Edges_keys"] = calc_dict["edges"]
        stats["Number_of_edges"] = len(calc_dict["edges"])
        stats["Number_of_node_types"] = len(calc_dict["types_freq"])
        return stats

    def describe(self) -> str:
        """
        Return a string describing the graph

        Returns:
        text based description of the graph
        """
        output = ""
        stats = self.statistics()
        for entry in stats:
            if "Number" in entry:
                string = str(entry) + "\t" + str(stats[entry])
                output = output + "\n" + string
        return output

    def get_upstream_branchpoint(self, node: str) -> list:
        """
        Return a list of branchpoints upstream of a node

        Parameters
        node: str
            key of the node

        Returns:
        list
        """
        if node in self.get_startpoints():
            return node

        upstream_node = self.nodes[node].input_nodes[0]

        if len(self.nodes[upstream_node].output_nodes) == 1:
            branchpoint = self.get_upstream_branchpoint(upstream_node)
            return branchpoint
        else:
            return upstream_node

    def newick(self, node=None, root=None) -> str:
        """
        Return a Newick representation of the graph.

        Parameters:
        node: Optional[str], default None
            Current node to process
        root: Optional[str], default None
            Root node of the tree

        Returns:
        str: Newick representation of the graph
        """
        # When the function is called initially without a node,
        # use the startpoint
        if node is None:
            startpoints = self.get_startpoints()
            if len(startpoints) > 1:
                raise ValueError("Graph has more than one startpoint")
            node = startpoints[0]
            root = node

        # Base case: If the current node is an endpoint,
        # return its label and branch length
        if node in self.get_endpoints():
            path_to_root = self.root_to_node_of_acyclic_node_path(node)
            start = self.nodes[node].position
            branch_length = start - self.nodes[path_to_root[-2]].position
            return f"{node}:{branch_length}"

        # Recursive case: Build the Newick string for
        # the children of the current node
        children = self.nodes[node].output_nodes

        # Iterate through the children and build the Newick string
        # for each of them
        child_newick_strings = []
        for child in children:
            child_newick = self.newick(child, root)
            child_newick_strings.append(child_newick)

        # Combine the Newick strings for the children with the current node
        newick = f"({','.join(child_newick_strings)}){node}"

        # Calculate branch length based on the difference between
        # node start and upstream node start
        if root != node:
            path_to_root = self.root_to_node_of_acyclic_node_path(node)
            start = self.nodes[node].position
            branch_length = start - self.nodes[path_to_root[-2]].position
            newick += f":{branch_length}"

        # If the current node is the root, append the final semicolon
        if node == root:
            newick += ";"

        return newick

    def newick_with_frames(self, node=None, root=None) -> str:
        """
        Return a Newick representation of the graph with frame information embedded in node labels.
        Preserves real branch lengths while adding frame information to node names.
        
        Parameters:
        node: Optional[str], default None
            Current node to process
        root: Optional[str], default None
            Root node of the tree
        Returns:
        str: Frame-enhanced Newick representation of the graph
        """
        # Create position to frame mapping from translons
        position_to_frame = {}
        translons = self.get_translons()
        
        for translon in translons:
            if len(translon) >= 2:
                start = translon[0]
                end = translon[1]
                frame = start % 3
                
                # Map all positions in this translon to this frame
                for pos in range(start, end + 1):
                    position_to_frame[pos] = frame
        
        def get_node_frame(node_key):
            """Get the reading frame for a node based on its position."""
            node_pos = self.nodes[node_key].position
            return position_to_frame.get(node_pos, -1)  # -1 for non-coding
        
        # When the function is called initially without a node,
        # use the startpoint
        if node is None:
            startpoints = self.get_startpoints()
            if len(startpoints) > 1:
                raise ValueError("Graph has more than one startpoint")
            node = startpoints[0]
            root = node
            
        # Base case: If the current node is an endpoint,
        # return its label with frame info and branch length
        if node in self.get_endpoints():
            path_to_root = self.root_to_node_of_acyclic_node_path(node)
            start = self.nodes[node].position
            branch_length = start - self.nodes[path_to_root[-2]].position
            frame = get_node_frame(node)
            frame_suffix = f"_F{frame}" if frame >= 0 else "_NC"  # NC for non-coding
            return f"{node}{frame_suffix}:{branch_length}"
            
        # Recursive case: Build the Newick string for
        # the children of the current node
        children = self.nodes[node].output_nodes
        
        # Iterate through the children and build the Newick string
        # for each of them
        child_newick_strings = []
        for child in children:
            child_newick = self.newick_with_frames(child, root)
            child_newick_strings.append(child_newick)
            
        # Get frame info for current node
        frame = get_node_frame(node)
        frame_suffix = f"_F{frame}" if frame >= 0 else "_NC"
        
        # Combine the Newick strings for the children with the current node
        newick = f"({','.join(child_newick_strings)}){node}{frame_suffix}"
        
        # Calculate branch length based on the difference between
        # node start and upstream node start
        if root != node:
            path_to_root = self.root_to_node_of_acyclic_node_path(node)
            start = self.nodes[node].position
            branch_length = start - self.nodes[path_to_root[-2]].position
            newick += f":{branch_length}"
            
        # If the current node is the root, append the final semicolon
        if node == root:
            newick += ";"
        return newick

    def prune(self, branch_position):
        """
        Prune the graph at the given position, removing all nodes and edges
        downstream of the position.

        Parameters:
        branch_position: int
            Position at which to prune the graph
        """
        affected_nodes = [
            node
            for node in self.nodes
            if self.nodes[node].position == branch_position
        ]

        for node in affected_nodes:
            for downstream_node in self.nodes[node].output_nodes:
                self.prune(self.nodes[downstream_node].position)

        for node in affected_nodes:
            self.remove_node(node)
