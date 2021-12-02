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


        self.node_count = 2 
        self.edge_count = 1 
        self.edges = edges
        self.nodes = nodes

    def load(self, locus_name="unnamed", locus_start=0, locus_stop=1000, nodes=None, edges=None):
        """
        create a RDG from the given paramaters 
        """
        if nodes:
            self.nodes = nodes
            self.node_count = len(nodes.keys())


        if edges:
            self.edges = edges 
            self.edge_count = len(edges.keys())
        

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
            1:Node(key=1, node_type="5_prime", coordinates=(0,1), edges_out=[1], nodes_out=[2]),
            2:Node(key=2, node_type="3_prime", coordinates=(1000 -1, 1000), edges_in=[2], nodes_in=[3]),
            3:Node(key=3, node_type="start", coordinates=(10,11), edges_in=[1], edges_out=[2,3], nodes_in=[1], nodes_out=[2,4]),
            4:Node(key=4, node_type="stop", coordinates=(100,101), edges_in=[3], edges_out=[4], nodes_in=[3], nodes_out=[5]),
            5:Node(key=5, node_type="3_prime", coordinates=(1000 -1, 1000), edges_in=[2], nodes_in=[3])
            }
        edges = {
            1:Edge(1, "untranslated", from_node=1, to_node=2, coordinates=(1, 1000 -1)),
            2:Edge(2, "untranslated", from_node=3, to_node=2, coordinates=(11, 1000 -1)),
            3:Edge(3, "translated", from_node=3, to_node=4, coordinates=(11, 100)),
            4:Edge(4, "untranslated", from_node=4, to_node=5, coordinates=(101, 1000 -1))
        }
        g = RDG()
        g = RDG.load(g, nodes=nodes, edges=edges)
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
        old_edge_coordinates = self.edges[edge_key].coordinates

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

        edge_key = self.get_new_edge_key()
        no_start = Edge(key=edge_key, edge_type="untranslated", from_node=start_node.key, to_node=edge.to_node, coordinates=(start_node.node_start, edge.coordinates[1]))

        self.add_edge(no_start, start_node.key, edge.to_node)


        edge_key = self.get_new_edge_key()
        three_prime = Edge(key=edge_key, edge_type="untranslated", from_node=stop_node.key, to_node=edge.to_node, coordinates=three_prime_edge_coords)

        terminal_node_key = self.get_new_node_key()
        terminal_node = Node(key=terminal_node_key, node_type="3_prime", coordinates=(edge.coordinates[1],edge.coordinates[1] + 1), edges_in=[], edges_out=[], nodes_in=[], nodes_out=[])
        self.add_node(terminal_node)
        self.add_edge(three_prime, stop_node.key, terminal_node.key)

        self.update_edge(edge.key, start_node.key, edge.to_node, (start_node.node_stop, self.nodes[edge.to_node].node_start))





    def add_open_reading_frame(self, start_codon_position, stop_codon_position):
        '''
        Handles all operations related to adding a new decision to the graph. Includes objects in graph and corrects all objects involved
        '''
        exisiting_edges = self.get_edges()
        clashing_edges = [] 
        node_key = self.get_new_node_key()
        start_node = Node(key=node_key, node_type="start_codon", coordinates=(start_codon_position, start_codon_position + 2), edges_in=[], edges_out=[], nodes_in=[], nodes_out=[])

        node_key = node_key + 1
        stop_node = Node(key=node_key, node_type="stop_codon", coordinates=(stop_codon_position, stop_codon_position + 2), edges_in=[], edges_out=[], nodes_in=[], nodes_out=[])

        for edge in exisiting_edges:
            if start_codon_position in range(self.edges[edge].coordinates[0], self.edges[edge].coordinates[1]):
                clashing_edges.append(edge)
        
        for edge in clashing_edges:
            self.insert_ORF(self.edges[edge], start_node, stop_node)
    
    def explore_node_ouputs(self, node, all_paths, current_path):
        '''
        explore the nodes directlty down path from input node 
        '''
        print(node, "current path", current_path, all_paths)
        if self.nodes[node].output_nodes:
            for output_node in self.nodes[node].output_nodes:
                current_path.append(output_node)
                if output_node == 2: print(node, self.nodes[node].output_nodes)

                if self.nodes[output_node].node_type == "3_prime":
                    print(current_path, all_paths)
                    all_paths.append(current_path)
                    current_path = []
                    return all_paths
                
                else:
                    # print("exploring: ", output_node, self.nodes[output_node].output_nodes)
                    self.explore_node_ouputs(output_node, all_paths, current_path)
        else:
            print(node, current_path)


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
        nodes = list(self.graph.keys())
        calc_dict = {
            "edges": [],
            "frames": [],
            "types": [],
            "heights": [],
            "frames_freq": {},
            "types_freq": {},
            "heights_freq": {},
        }

        for i in nodes:
            for j in self.graph[i].connected:
                if (j, i) not in calc_dict["edges"]:
                    calc_dict["edges"].append((i, j))

            calc_dict["frames"].append(self.graph[i].frame)
            calc_dict["types"].append(self.graph[i].node_type)
            calc_dict["heights"].append(self.graph[i].height)

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

        for item in calc_dict["heights"]:
            if item in calc_dict["heights_freq"]:
                calc_dict["heights_freq"][item] += 1
            else:
                calc_dict["heights_freq"][item] = 1

        for item in calc_dict["types_freq"]:
            key = "Number_of_" + item
            stats[key] = calc_dict["types_freq"][item]

        for item in calc_dict["frames_freq"]:
            key = "Number_of_nodes_in_frame_" + str(item)
            stats[key] = calc_dict["frames_freq"][item]

        for item in calc_dict["heights_freq"]:
            key = "Number_of_nodes_with_height_" + str(item)
            stats[key] = calc_dict["heights_freq"][item]

        stats["Node_keys"] = nodes
        stats["Number_of_nodes"] = self.node_count
        stats["Edges_keys"] = calc_dict["edges"]
        stats["Number_of_edges"] = len(calc_dict["edges"])
        stats["Number_of_node_types"] = len(calc_dict["types_freq"])

        return stats

    def describe(self):
        """
        print in the terminal a description of the graph
        """

        stats = self.statistics()
        for entry in stats:
            print(str(entry) + "\t" + str(stats[entry]))





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

    g = { "a" : {"out":[], "in":[]},
          "b" : {"out":[], "in":[]},
          "c" : {"out":[], "in":[]},
          "d" : {"out":[], "in":[]},
          "e" : {"out":[], "in":[]},
          "f" : {"out":[], "in":[]}
        }

    dg = RDG()
    dg = dg.load_example()
    dg.remove_edge(1)
    # dg.add_open_reading_frame(30, 90)

    # dg.add_open_reading_frame(150, 171)
    # dg.print_paths()


    # for i in dg.nodes:
        # print(dg.nodes[i].key, dg.nodes[i].node_start, dg.nodes[i].node_stop, dg.nodes[i].node_type, dg.nodes[i].output_nodes)
    # print ()

    # for i in dg.edges:
    #     print(dg.edges[i].key, dg.edges[i].coordinates, dg.edges[i].edge_type)

