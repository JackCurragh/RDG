
class Node(object):
    def __init__(self, key, node_type, position, edges_in=[], edges_out=[], nodes_in=[], nodes_out=[]):
        self.key = key

        valid_node_types = ["5_prime", "3_prime", "start", "stop", "frameshift", "readthrough_stop"]
        if node_type not in valid_node_types:
            raise ValueError(f"Node type ({node_type}) not valid. Valid types are {valid_node_types}")
        
        self.node_type = node_type
        self.input_edges = edges_in
        self.output_edges = edges_out
        self.input_nodes = nodes_in 
        self.output_nodes = nodes_out 
        self.node_start = position
        self.frame = self.node_start % 3

    def get_node_key(self):
        return self.key

    def get_node_type(self):
        return self.node_type


class Edge(object):
    def __init__(self, key, edge_type, from_node, to_node, coordinates):
        self.key = key
        self.edge_type = edge_type
        self.from_node = from_node
        self.to_node = to_node
        self.coordinates = coordinates

    def edge_key(self):
        return self.key

    def edge_type(self):
        return self.edge_type

    def from_node(self):
        return self.from_node

    def to_node(self):
        return self.to_node 
    
    def coordinates(self):
        return self.coordinates

    def frame(self):
        return self.coordinates[0] % 3 



class RDG(object):
    def __init__(self, locus_stop: int=1000, name: str=None):
        """
        Initialise a fresh RDG object 

        Parameters:
        -----------
        locus_stop: int length of the locus in bp
        name: str name of the locus

        """
        self.locus = name
        self.locus_start = 0
        self.locus_stop = locus_stop

        transcript_length = locus_stop - self.locus_start

        self.nodes = {
                1:Node(key=1, node_type="5_prime", position=0, edges_out=[1], nodes_out=[2]),
                2:Node(key=2, node_type="3_prime", position=transcript_length - 1, edges_in=[1], nodes_in=[1])
            }

        self.edges = {1: Edge(1, "untranslated", from_node=1, to_node=2, coordinates=(1, transcript_length -1))}


    def Load(self, locus_name: str="unnamed", locus_start: int=0, locus_stop: int=1000, nodes: dict=None, edges: dict=None):
        """
        Create a RDG from the given paramaters 

        Parameters:
        -----------
        locus_name: str name of the locus
        locus_start: int start of the locus in bp
        locus_stop: int stop of the locus in bp
        nodes: dict of nodes. Keys are node keys, values are Node objects
        edges: dict of edges. Keys are edge keys, values are Edge objects

        Returns:
        --------
        RDG object
        """

        if not nodes:
            raise ValueError("nodes and edges must be provided to load a graph. If you want to create a new graph, use RDG()")
        
        if locus_stop < locus_start:
            raise ValueError("locus_stop must be greater than locus_start")

        self.nodes = nodes
        self.edges = edges 
        self.locus = locus_name
        self.locus_start = locus_start
        self.locus_stop = locus_stop
        return self


    def load_example(self):
        '''
        Load a basic graph with one ORF 

        Returns:
        --------
        RDG object
        with:
            locus length = 1000
            orf coordinates = 10, 100
        '''
        nodes = {
            1:Node(key=1, node_type="5_prime", position=0, edges_out=[1], nodes_out=[3]),
            2:Node(key=2, node_type="3_prime", position=1000 -1, edges_in=[2], nodes_in=[3]),
            3:Node(key=3, node_type="start", position=10, edges_in=[1], edges_out=[2,3], nodes_in=[1], nodes_out=[2,4]),
            4:Node(key=4, node_type="stop", position=100, edges_in=[3], edges_out=[4], nodes_in=[3], nodes_out=[5]),
            5:Node(key=5, node_type="3_prime", position=1000 -1, edges_in=[4], nodes_in=[4])
            }
        edges = {
            1:Edge(1, "untranslated", from_node=1, to_node=3, coordinates=(1, 10 -1)),
            2:Edge(2, "untranslated", from_node=3, to_node=2, coordinates=(11, 1000 -1)),
            3:Edge(3, "translated", from_node=3, to_node=4, coordinates=(11, 100)),
            4:Edge(4, "untranslated", from_node=4, to_node=5, coordinates=(101, 1000 -1))
        }
        g = RDG()
        g = RDG.Load(g, nodes=nodes, edges=edges)
        return g


    def get_nodes(self) -> list:
        """
        Return a list of the nodes from the graph

        Returns:
        --------
        list of Node ids
        """
        return list(self.nodes.keys())


    def get_edges(self) -> list:
        '''
        Return a list of keys of edges in the graph

        Returns:
        --------
        list of Edge ids
        '''
        return list(self.edges.keys())


    def get_edges_from_to(self) -> dict:
        '''
        Return a dict of edges with keys of the form (from_node_id, to_node_id)

        Returns:
        --------
        dict Keys are tuples of the form (from_node_id, to_node_id) Values are edge ids
        '''
        return {(self.edges[edge].from_node, self.edges[edge].to_node): edge for edge in self.edges.keys()}


    def get_new_node_key(self) -> int:
        """
        Return the lowest available node key

        Returns:
        --------
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
        --------
        int: new edge key that is not already used
        """
        keys = self.edges.keys()
        if len(keys) > 0:
            return max(keys) + 1
        else:
            return 1


    def get_key_from_position(self, position, node_type) -> int:
        '''
        Return the node key for the node of specified type at the stated position 

        Parameters:
        -----------
        position: int position in the locus
        node_type: str type of node to look for

        Returns:
        --------
        int: node key

        '''

        valid_nodes = {node: self.nodes[node] for node in self.nodes if self.nodes[node].node_type == node_type}

        if valid_nodes == {}:
            raise ValueError(f"There are no nodes of type '{node_type}' in the graph")
        else:
            for node in valid_nodes:
                if valid_nodes[node].node_start == position:
                    return node
            raise ValueError(f"There is no node of type '{node_type}' at position {position}")




    def remove_edge(self, edge_key):
        '''
        Remove the edge with the given key from the graph also removing references from adjacent nodes 

        Does not remove nodes that are no longer connected to any other nodes

        Parameters:
        -----------
        edge_key: int key of the edge to remove

        '''
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


    def add_edge(self, edge, from_node_key, to_node_key):
        '''
        Introduce edges and references with the a edge 

        Parameters:
        -----------
        edge: Edge object to add
        from_node_key: int key of the node the edge is coming from
        to_node_key: int key of the node the edge is going to

        '''
        self.nodes[from_node_key].output_edges.append(edge.key)
        self.nodes[from_node_key].output_nodes.append(to_node_key)

        self.nodes[to_node_key].input_edges.append(edge.key)
        self.nodes[to_node_key].input_nodes.append(from_node_key)

        self.edges[edge.key] = edge 
        self.edges[edge.key].frame = self.nodes[from_node_key].node_start%3



    def add_node(self, node):
        '''
        Introduce nodes and references with the a node. 

        The associated edges must be added separately

        Parameters:
        -----------
        node: Node object to add

        '''
        self.nodes[node.key] = node


    def update_edge(self, edge_key, new_from_node, new_to_node, new_coordinates):
        '''
        Update the edge with the given key with the new from and to nodes and coordinates

        Parameters:
        -----------
        edge_key: int key of the edge to update
        new_from_node: int key of the new from node
        new_to_node: int key of the new to node
        new_coordinates: tuple of the form (start, end) of the new coordinates

        notes:
        when inserting an orf into an untranslated edge we update the untranslated edge to now start at the new start codon instead

        '''
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

        if old_to_node != new_to_node:
            if edge_key in self.nodes[old_to_node].input_edges:
                self.nodes[old_to_node].input_edges.remove(edge_key)

            if edge_key not in self.nodes[new_to_node].input_edges:
                self.nodes[new_to_node].input_edges.append(edge_key)

            if old_from_node in self.nodes[old_to_node].input_nodes:
                self.nodes[old_to_node].input_nodes.remove(old_from_node)

            if new_from_node in self.nodes[new_to_node].input_nodes:
                self.nodes[new_to_node].input_nodes.remove(new_from_node)

        self.edges[edge_key].coordinates = new_coordinates
        self.edges[edge_key].from_node = new_from_node


    def insert_ORF(self, edge: Edge, start_node: Node, stop_node: Node):
        '''
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

        '''
        five_prime_edge_coords = (edge.coordinates[0], start_node.node_start - 1)
        coding_edge_coords = (start_node.node_start, stop_node.node_start)
        three_prime_edge_coords = (stop_node.node_start + 1, edge.coordinates[1])

        edge_key = self.get_new_edge_key()
        five_prime = Edge(key=edge_key, edge_type="untranslated", from_node=edge.from_node, to_node=start_node.key, coordinates=five_prime_edge_coords)

        self.add_node(start_node)
        self.add_edge(five_prime, edge.from_node, start_node.key)

        edge_key = self.get_new_edge_key()
        coding = Edge(key=edge_key, edge_type="translated", from_node=start_node.key, to_node=stop_node.key, coordinates=coding_edge_coords)

        self.add_node(stop_node)
        self.add_edge(coding, start_node.key, stop_node.key)


        terminal_node_key = self.get_new_node_key()
        terminal_node = Node(key=terminal_node_key, node_type="3_prime", position=edge.coordinates[1], edges_in=[], edges_out=[], nodes_in=[], nodes_out=[])

        edge_key = self.get_new_edge_key()
        three_prime = Edge(key=edge_key, edge_type="untranslated", from_node=stop_node.key, to_node=terminal_node.key, coordinates=three_prime_edge_coords)

        self.add_node(terminal_node)
        self.add_edge(three_prime, stop_node.key, terminal_node.key)

        self.update_edge(edge.key, start_node.key, edge.to_node, (start_node.node_start, self.nodes[edge.to_node].node_start))


    def is_input_edge_translated(self, node: int) -> bool:
        '''
        Return boolean as to whether an edge entering this node is translated 

        Parameters:
        -----------
        node: int key of the node to check

        Returns:
        --------
        boolean: True if an edge entering this node is translated, False otherwise
        '''
        input_edges = self.nodes[node].input_edges
        boolean = False
        for edge in input_edges:
            if self.edges[edge].edge_type == 'translated':
                boolean = True
        
        return boolean


    def root_to_node_of_acyclic_node_path(self, node) -> list:
        '''
        Return path from root to given node as a list of node keys

        Parameters:
        -----------
        node: int key of the node to check

        Returns:
        --------
        path: list of node keys from root to given node
        '''
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


    def root_to_node_of_acyclic_edge_path(self, node):
        '''
        Return path from root to given node as a list of edge keys

        Parameters:
        -----------
        node: int key of the node to check

        Returns:
        --------
        path: list of edge keys from root to given node

        '''
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
            

    def check_translation_upstream(self, from_node, upstream_limit=1):
        '''
        look upstream of an edges from node and see if any edges are of type 'translated' 
        used in reinit functionality. Upstream limit refers to the number of ORFs allowed upstream. ie. can the be multiple reinitiation events or just one etc
        '''
        edge_path_to_node_from_root = self.root_to_node_of_acyclic_edge_path(from_node)
        number_of_translated_regions = 0
        for edge in edge_path_to_node_from_root:
            if self.edges[edge].edge_type == 'translated':
                number_of_translated_regions += 1

        if number_of_translated_regions <= upstream_limit:
            return False
        return True


    def add_open_reading_frame(self, start_codon_position, stop_codon_position, reinitiation=False, upstream_limit=0):
        '''
        Handles all operations related to adding a new decision to the graph. Includes objects in graph and corrects all objects involved
        '''
        exisiting_edges = self.get_edges()
        clashing_edges = [] 
        for edge in exisiting_edges:
            if start_codon_position in range(self.edges[edge].coordinates[0], self.edges[edge].coordinates[1]) and self.edges[edge].edge_type != "translated":
                upstream_node = self.edges[edge].from_node
                clashing_edges.append((edge, upstream_node))
                
        
        for edge, upstream_node in clashing_edges:
            node_key = self.get_new_node_key()
            start_node = Node(key=node_key, node_type="start", position=start_codon_position, edges_in=[], edges_out=[], nodes_in=[], nodes_out=[])

            node_key = node_key + 1
            stop_node = Node(key=node_key, node_type="stop", position=stop_codon_position, edges_in=[], edges_out=[], nodes_in=[], nodes_out=[])
            if reinitiation or not self.check_translation_upstream(upstream_node, upstream_limit=upstream_limit):
                self.insert_ORF(self.edges[edge], start_node, stop_node)

    
    def add_stop_codon_readthrough(self, readthrough_codon_position, next_stop_codon_position):
        '''
        Handles all operations related to adding a new decision to the graph at a stop codon. Includes objects in graph and corrects all objects involved

        '''
        readthrough_codon_key = self.get_key_from_position(readthrough_codon_position, "stop")
        three_prime_terminal_key = self.nodes[readthrough_codon_key].output_nodes[0]

        node_key = self.get_new_node_key()
        terminal_node_key = node_key + 1

        edge_key = self.get_new_edge_key()
        edge_key2 = edge_key + 1

        new_stop_node = Node(key=node_key, node_type="stop", position=next_stop_codon_position, edges_in=[edge_key], edges_out=[edge_key2], nodes_in=[readthrough_codon_key], nodes_out=[terminal_node_key])
        self.add_node(new_stop_node)

        coding = Edge(key=edge_key, edge_type="translated", from_node=readthrough_codon_key, to_node=new_stop_node.key, coordinates=(self.nodes[readthrough_codon_key].node_start, next_stop_codon_position))
        self.add_edge(coding, readthrough_codon_key, new_stop_node.key)

        self.nodes[readthrough_codon_key].output_edges.append(coding.key)
        self.nodes[readthrough_codon_key].output_nodes.append(new_stop_node.key)
        self.nodes[readthrough_codon_key].node_type = "readthrough_stop"
        
        three_prime = Edge(key=edge_key2, edge_type="untranslated", from_node=new_stop_node.key, to_node=terminal_node_key, coordinates=(new_stop_node.node_start, self.nodes[three_prime_terminal_key].node_start))

        terminal_node_key = self.get_new_node_key()
        terminal_node = Node(key=terminal_node_key, node_type="3_prime", position=self.nodes[three_prime_terminal_key].node_start, edges_in=[edge_key2], edges_out=[], nodes_in=[new_stop_node.key], nodes_out=[])
        self.add_node(terminal_node)

        self.add_edge(three_prime, new_stop_node.key, terminal_node_key)
        self.add_edge(three_prime, readthrough_codon_key, three_prime_terminal_key)


    def add_frameshift(self, fs_position, next_stop_codon_position, shift):
        '''
        add a frameshifting event 
        '''
        exisiting_edges = self.get_edges()
        clashing_edges = [] 
        for edge in exisiting_edges:
            if fs_position in range(self.edges[edge].coordinates[0], self.edges[edge].coordinates[1]) and self.edges[edge].edge_type == "translated":
                upstream_node = self.edges[edge].from_node
                clashing_edges.append((edge, upstream_node))  

        for edge, upstream_node in clashing_edges:
            # edge and upstrream are keys
            shift_node_key = self.get_new_node_key()
            new_stop_node_key = shift_node_key + 1
            old_stop_node_key = self.edges[edge].to_node
            new_terminal_key = new_stop_node_key + 1

            old_stop_edge_key = self.get_new_edge_key()
            new_stop_edge_key = old_stop_edge_key + 1
            new_three_prime_edge_key = new_stop_edge_key + 1

            shift_node = Node(key=shift_node_key, 
            node_type="frameshift", 
            coordinates=fs_position, 
            edges_in=[edge], 
            edges_out=[old_stop_edge_key, new_stop_edge_key], 
            nodes_in=[upstream_node], 
            nodes_out=[old_stop_node_key, new_stop_node_key])
            self.add_node(shift_node)

            self.nodes[upstream_node].output_nodes.remove(old_stop_node_key)
            self.nodes[upstream_node].output_nodes.append(shift_node.key)

            self.edges[edge].to_node = shift_node_key

            self.nodes[old_stop_node_key].input_edges.remove(edge)
            self.nodes[old_stop_node_key].input_edges.append(old_stop_edge_key)
            self.nodes[old_stop_node_key].input_nodes.remove(upstream_node)
            self.nodes[old_stop_node_key].input_nodes.append(shift_node_key)

            
            new_stop_node = Node(key=new_stop_node_key, 
            node_type="stop", 
            coordinates=next_stop_codon_position, 
            edges_in=[new_stop_edge_key], 
            edges_out=[new_three_prime_edge_key], 
            nodes_in=[shift_node_key], 
            nodes_out=[new_terminal_key])
            self.add_node(new_stop_node)
            
            new_terminal_node = Node(key=new_terminal_key, 
            node_type="3_prime",
            coordinates=self.locus_stop, 
            edges_in=[new_three_prime_edge_key], 
            edges_out=[], 
            nodes_in=[new_stop_node_key], 
            nodes_out=[])
            self.add_node(new_terminal_node)

            old_stop_edge = Edge(key=old_stop_edge_key,
            edge_type="translated",
            from_node=shift_node_key,
            to_node=old_stop_node_key,
            coordinates=(shift_node.node_start, self.nodes[old_stop_node_key].node_start))
            self.add_edge(old_stop_edge, from_node_key=shift_node_key, to_node_key=old_stop_node_key)
            self.edges[old_stop_edge_key].frame = self.edges[edge].frame

            new_stop_edge = Edge(key=new_stop_edge_key,
            edge_type="translated",
            from_node=shift_node_key,
            to_node=new_stop_node_key,
            coordinates=(shift_node.node_start, self.nodes[new_stop_node_key].node_start))
            self.add_edge(new_stop_edge, from_node_key=shift_node_key, to_node_key=new_stop_node_key)


            new_three_prime_edge = Edge(key=new_three_prime_edge_key,
            edge_type="untranslated",
            from_node=new_stop_node_key,
            to_node=new_terminal_key,
            coordinates=(new_stop_node.node_start, new_terminal_node.node_start))
            self.add_edge(new_three_prime_edge, from_node_key=new_stop_node_key, to_node_key=new_terminal_key)



    def get_branch_points(self):
        '''
        return a list of nodes that are branch points
        '''
        branch_points = []
        for node in self.nodes:
            if len(self.nodes[node].output_edges) > 1:
                branch_points.append(node)
        return branch_points


    def get_endpoints(self):
        '''
        return list of terminal node keys 
        '''
        endpoints = [] 
        for node in self.nodes:
            if len(self.nodes[node].output_edges) == 0 and self.nodes[node].node_type == "3_prime":
                endpoints.append(node)
        return endpoints


    def get_startpoints(self):
        '''
        return list of start node keys 
        '''
        startpoints = [] 
        for node in self.nodes:
            if len(self.nodes[node].input_edges) == 0 and self.nodes[node].node_type == "5_prime":
                startpoints.append(node)
        return startpoints
    

    def get_start_nodes(self):
        '''
        return translation start keys
        '''
        translation_starts = [] 
        for node in self.nodes:
            if len(self.nodes[node].output_edges) == 2 and self.nodes[node].node_type == "start":
                translation_starts.append(node)
        return translation_starts


    def get_stop_nodes(self):
        '''
        return translation start keys
        '''
        translation_stops = [] 
        for node in self.nodes:
            if "stop" in  self.nodes[node].node_type:
                translation_stops.append(node)
        return translation_stops

    
    def get_orfs(self):
        '''
        return a list of coordinates (start, stop) describing the ORFs in the graph
        ''' 
        translation_starts = self.get_start_nodes()

        orfs = []
        for start in translation_starts:
            downstream_nodes = self.nodes[start].output_nodes

            for candidate_stop in downstream_nodes:
                if self.nodes[candidate_stop].node_type == "stop":
                    orfs.append((self.nodes[start].node_start, self.nodes[candidate_stop].node_start))

                elif self.nodes[candidate_stop].node_type == "frameshift":
                    # orfs.append((self.nodes[start].node_start, self.nodes[candidate_stop].node_start))
                    fs_downstream_nodes = self.nodes[candidate_stop].output_nodes
                    for new_candidate_stop in fs_downstream_nodes:
                        if self.nodes[new_candidate_stop].node_type == "stop":
                            orf = (self.nodes[start].node_start, self.nodes[candidate_stop].node_start, self.nodes[new_candidate_stop].node_start)
                            if orf not in orfs:
                                orfs.append(orf)

                elif self.nodes[candidate_stop].node_type == "readthrough_stop":
                    orfs.append((self.nodes[start].node_start, self.nodes[candidate_stop].node_start))
                    readthrough_downstream_nodes = self.nodes[candidate_stop].output_nodes

                    for new_candidate_stop in readthrough_downstream_nodes:
                        if self.nodes[new_candidate_stop].node_type == "stop":
                            orf = (self.nodes[start].node_start, self.nodes[new_candidate_stop].node_start)
                            if orf not in orfs:
                                orfs.append(orf)


        
        return orfs


    def get_frameshifts(self):
        '''
        return a list of keys for all frameshift nodes
        '''
        frameshifts = [] 
        for node in self.nodes:
            if "frameshift" in  self.nodes[node].node_type:
                frameshifts.append(node)
        return frameshifts


    def print_paths(self):
        '''
        read out the paths from the graph object 
        '''
        root_nodes = [] 
        for node in self.nodes:
            if self.nodes[node].node_type == "5_prime":
                root_nodes.append(node)

        all_paths = [] 
        for root in root_nodes:
            path_from_root_nodes = [root] 
            if self.nodes[root].output_nodes:
                for node in self.nodes[root].output_nodes:
                    self.explore_node_ouputs(node, all_paths, path_from_root_nodes)

            else:
                all_paths.append(path_from_root_nodes)
        print("all paths", all_paths)


    def statistics(self):
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
            key = "Number_of_nodes_in_frame_" + str(item)
            stats[key] = calc_dict["frames_freq"][item]


        stats["Node_keys"] = nodes
        stats["Number_of_nodes"] = len(self.nodes.keys())
        stats["Edges_keys"] = calc_dict["edges"]
        stats["Number_of_edges"] = len(calc_dict["edges"])
        stats["Number_of_node_types"] = len(calc_dict["types_freq"])
        return stats


    def describe(self):
        """
        print in the terminal a description of the graph
        """
        output = ""
        stats = self.statistics()
        for entry in stats:
            if "Number" in entry: 
                string = str(entry) + "\t" + str(stats[entry])
                output = output  + "\n" + string 
        return output









# if __name__ == "__main__":

#     dg = RDG()
#     dg = dg.load_example()
#     dg.add_open_reading_frame(30, 90)

#     dg.add_open_reading_frame(150, 171)
#     dg.add_stop_codon_readthrough(90, 120)
#     print(dg.get_orfs())
