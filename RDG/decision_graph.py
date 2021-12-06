class RDG(object):
    def __init__(self, locus_stop=1000, nodes=None, edges=None):
        """
        initialise a fresh empty RDG 
        """
        self.locus = ""
        self.locus_start = 0
        self.locus_stop = locus_stop

        if nodes: 
            raise ValueError("A value for 'nodes' was provided. If you want to specify from an existing structure use RDG.load()")
        if edges: 
            raise ValueError("A value for 'edges' was provided. If you want to specify from an existing structure use RDG.load()")

        transcript_length = locus_stop - self.locus_start

        nodes = {
                1:Node(key=1, node_type="5_prime", coordinates=(0,1), edges_out=[1], nodes_out=[2]),
                2:Node(key=2, node_type="3_prime", coordinates=(transcript_length -1, transcript_length), edges_in=[1], nodes_in=[1])
            }

        edges = {1: Edge(1, "untranslated", from_node=1, to_node=2, coordinates=(1, transcript_length -1))}


        self.edges = edges
        self.nodes = nodes

    def Load(self, locus_name="unnamed", locus_start=0, locus_stop=1000, nodes=None, edges=None):
        """
        create a RDG from the given paramaters 
        """
        if nodes:
            self.nodes = nodes

        if edges:
            self.edges = edges 
        
        self.locus = locus_name
        self.locus_start = locus_start
        self.locus_stop = locus_stop
        return self

    def load_example(self):
        '''
        load a basic graph with one ORF 
        locus length = 1000
        orf coordinates = 10, 100
        '''
        nodes = {
            1:Node(key=1, node_type="5_prime", coordinates=(0,1), edges_out=[1], nodes_out=[3]),
            2:Node(key=2, node_type="3_prime", coordinates=(1000 -1, 1000), edges_in=[2], nodes_in=[3]),
            3:Node(key=3, node_type="start", coordinates=(10,11), edges_in=[1], edges_out=[2,3], nodes_in=[1], nodes_out=[2,4]),
            4:Node(key=4, node_type="stop", coordinates=(100,101), edges_in=[3], edges_out=[4], nodes_in=[3], nodes_out=[5]),
            5:Node(key=5, node_type="3_prime", coordinates=(1000 -1, 1000), edges_in=[4], nodes_in=[4])
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


    def vertices(self):
        """
        return a list of the vertices/nodes from the graph
        """
        return list(self.nodes.keys())


    def get_edges(self):
        '''
        return a list of keys of edges in the graph
        '''
        return list(self.edges.keys())


    def get_edges_from_to(self):
        '''
        return a list of edges (from, to) in the graph
        '''
        from_to = {} 
        for edge in self.edges.keys():
            from_to[(self.edges[edge].from_node, self.edges[edge].to_node)] = edge

        return from_to


    def get_new_node_key(self):
        """
        check the nodes to find the lowest available key
        """
        keys = self.nodes.keys()
        if len(keys) > 0:
            return max(keys) + 1
        else:
            return 1


    def get_new_edge_key(self):
        """
        check the edges to find the lowest available key
        """
        keys = self.edges.keys()
        if len(keys) > 0:
            return max(keys) + 1
        else:
            return 1

    def remove_edge(self, edge_key):
        '''
        remove the edge with the given key from the graph Also removing references from adjacent nodes 
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
        reintroduce edges and references with the new edge 
        '''
        self.nodes[from_node_key].output_edges.append(edge.key)
        self.nodes[from_node_key].output_nodes.append(to_node_key)

        self.nodes[to_node_key].input_edges.append(edge.key)
        self.nodes[to_node_key].input_nodes.append(from_node_key)

        self.edges[edge.key] = edge 


    def add_node(self, node):
        self.nodes[node.key] = node


    def update_edge(self, edge_key, new_from_node, new_to_node, new_coordinates):
        '''
        functionality for changing 
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


    def insert_ORF(self, edge, start_node, stop_node):
        '''
        handle inserting node into just one edge including down path 
        '''
        five_prime_edge_coords = (edge.coordinates[0], start_node.node_start - 1)
        coding_edge_coords = (start_node.node_start, stop_node.node_stop)
        three_prime_edge_coords = (stop_node.node_stop + 1, edge.coordinates[1])

        edge_key = self.get_new_edge_key()
        five_prime = Edge(key=edge_key, edge_type="untranslated", from_node=edge.from_node, to_node=start_node.key, coordinates=five_prime_edge_coords)

        self.add_node(start_node)
        self.add_edge(five_prime, edge.from_node, start_node.key)

        edge_key = self.get_new_edge_key()
        coding = Edge(key=edge_key, edge_type="translated", from_node=start_node.key, to_node=stop_node.key, coordinates=coding_edge_coords)

        self.add_node(stop_node)
        self.add_edge(coding, start_node.key, stop_node.key)


        terminal_node_key = self.get_new_node_key()
        terminal_node = Node(key=terminal_node_key, node_type="3_prime", coordinates=(edge.coordinates[1],edge.coordinates[1] + 1), edges_in=[], edges_out=[], nodes_in=[], nodes_out=[])

        edge_key = self.get_new_edge_key()
        three_prime = Edge(key=edge_key, edge_type="untranslated", from_node=stop_node.key, to_node=terminal_node.key, coordinates=three_prime_edge_coords)

        self.add_node(terminal_node)
        self.add_edge(three_prime, stop_node.key, terminal_node.key)

        self.update_edge(edge.key, start_node.key, edge.to_node, (start_node.node_stop, self.nodes[edge.to_node].node_start))


    def is_input_edge_translated(self, node):
        '''
        return boolean as to whether an edge entering this node is translated 
        '''
        input_edges = self.nodes[node].input_edges
        boolean = False
        for edge in input_edges:
            if self.edges[edge].edge_type == 'translated':
                boolean = True
        
        return boolean


    def root_to_node_of_acyclic_node_path(self, node):
        '''
        return path from root to given node as a list of node keys
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
        return(path[::-1])


    def root_to_node_of_acyclic_edge_path(self, node):
        '''
        return path from root to given node as a list of edge keys
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
            

    def check_translation_upstream(self, from_node):
        '''
        look upstream of an edges from node and see if any edges are of type 'translated' 
        used in reinit functionality 
        '''
        edge_path_to_node_from_root = self.root_to_node_of_acyclic_edge_path(from_node)

        for edge in edge_path_to_node_from_root:
            if self.edges[edge].edge_type == 'translated':
                return True
        return False

    def add_open_reading_frame(self, start_codon_position, stop_codon_position, reinitiation=False):
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
            start_node = Node(key=node_key, node_type="start", coordinates=(start_codon_position, start_codon_position + 2), edges_in=[], edges_out=[], nodes_in=[], nodes_out=[])

            node_key = node_key + 1
            stop_node = Node(key=node_key, node_type="stop", coordinates=(stop_codon_position, stop_codon_position + 2), edges_in=[], edges_out=[], nodes_in=[], nodes_out=[])
            if reinitiation or not self.check_translation_upstream(upstream_node):
                self.insert_ORF(self.edges[edge], start_node, stop_node)
    

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
            if self.nodes[node].node_type == "stop":
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
        
        return orfs






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







class Node(object):
    def __init__(self, key, node_type, coordinates, edges_in=[], edges_out=[], nodes_in=[], nodes_out=[]):
        self.key = key
        self.node_type = node_type
        self.input_edges = edges_in
        self.output_edges = edges_out
        self.input_nodes = nodes_in 
        self.output_nodes = nodes_out 
        self.node_start = coordinates[0]
        self.node_stop = coordinates[1]
        self.frame = self.node_start % 3


    def add_neighbour(self, neighbour, weight=0):
        self.connected[neighbour] = weight

    def get_connections(self):
        return self.connected.keys()

    def node_key(self):
        return self.key

    def node_type(self):
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



if __name__ == "__main__":

    dg = RDG()
    dg = dg.load_example()
    dg.add_open_reading_frame(30, 90)

    dg.add_open_reading_frame(150, 171)
    print(dg.get_orfs())
