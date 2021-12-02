import matplotlib.pyplot as plt
import networkx as nx
from RDG import RDG, Node, Edge



def plot(graph):
    G = nx.Graph()
    edges = graph.get_edges_from_to()
    print(edges)
    G.add_edges_from(edges)
    nx.draw_networkx(G)
    plt.show()

    for i in graph.nodes:
        print(i)

def plot_directed(graph):
    G = nx.DiGraph()
    edges = graph.get_edges_from_to()

    G.add_edges_from(edges)
    nx.draw_networkx(G)
    plt.show()


if __name__ == "__main__":
    dg = RDG()
    dg = dg.load_example()
    node_key = dg.get_new_node_key()

    start_codon_position = 15
    stop_codon_position = 25
    start_node = Node(key=node_key, node_type="start_codon", coordinates=(start_codon_position, start_codon_position + 2), edges_in=[], edges_out=[], nodes_in=[], nodes_out=[])

    node_key = node_key + 1
    stop_node = Node(key=node_key, node_type="stop_codon", coordinates=(stop_codon_position, stop_codon_position + 2), edges_in=[], edges_out=[], nodes_in=[], nodes_out=[])
    
    dg.insert_ORF(dg.edges[2], start_node, stop_node)
    plot_directed(dg)
