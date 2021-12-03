import matplotlib.pyplot as plt
import networkx as nx
from RDG import RDG, Node, Edge



def plot(graph):
    '''
    Plot a riboosome decision graph as a basic undirected graph with spring layout
    '''
    G = nx.Graph()
    edges = graph.get_edges_from_to()
    print(edges)
    G.add_edges_from(edges)
    nx.draw_networkx(G)
    plt.show()

    for i in graph.nodes:
        print(i)


def position_graphs_nodes(graph):
    '''
    determine the positioning of rdg nodes from graph structure. 
    This is a hacky method that determines Y axis by node key rather than calculating an optimum based on overlaps
    '''

    node_x_positions = [(graph.nodes[node].key, graph.nodes[node].node_start) for node in graph.nodes]
    pos = {}
    for i in node_x_positions:
        pos[i[0]] = (i[1], i[0])

    return pos

def plot_directed(graph):
    G = nx.DiGraph()
    edges = graph.get_edges_from_to()
    pos = position_graphs_nodes(graph)

    endpoints = graph.get_endpoints()
    startpoints = graph.get_startpoints()
    translation_starts = graph.get_start_nodes()
    translation_stops = graph.get_stop_nodes()

    G.add_edges_from(edges)
    node_colors = []
    for node in G.nodes():
        if node in startpoints:
            node_colors.append((0, 0, 1))
        elif node in endpoints:
            node_colors.append((0.5,0,0.5))
        elif node in translation_starts:
            node_colors.append((0,1,0))
        elif node in translation_stops:
            node_colors.append((1,0,0))
        else:
            node_colors.append((0,0,0))
    nx.draw_networkx(G, pos=pos, node_shape='s', node_color=node_colors)

    plt.show()


if __name__ == "__main__":
    dg = RDG()
    dg = dg.load_example()
    node_key = dg.get_new_node_key()

    start_codon_position = 150
    stop_codon_position = 850
    start_node = Node(key=node_key, node_type="start_codon", coordinates=(start_codon_position, start_codon_position + 2), edges_in=[], edges_out=[], nodes_in=[], nodes_out=[])

    node_key = node_key + 1
    stop_node = Node(key=node_key, node_type="stop_codon", coordinates=(stop_codon_position, stop_codon_position + 2), edges_in=[], edges_out=[], nodes_in=[], nodes_out=[])
    
    # dg.insert_ORF(dg.edges[2], start_node, stop_node)
    dg.add_open_reading_frame(150,850)
    plot_directed(dg)
