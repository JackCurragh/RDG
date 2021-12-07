import matplotlib.pyplot as plt
import networkx as nx
from RDG import RDG, Node, Edge



def position_graphs_nodes(graph):
    '''
    determine the positioning of rdg nodes from graph structure. 
    This is a hacky method that determines Y axis by node key rather than calculating an optimum based on overlaps
    '''

    node_x_positions = [(graph.nodes[node].key, graph.nodes[node].node_start) for node in graph.nodes]
    pos = {}
    for i in node_x_positions:
        pos[i[0]] = (i[1], i[0])
    
    highest_node = max(list(pos.keys()))
    pos[2] = (pos[2][0], pos[highest_node][1] + 1)

    return pos


def plot(graph, color_dict=None):
    G = nx.DiGraph()
    edges = graph.get_edges_from_to()

    pos = position_graphs_nodes(graph)

    endpoints = graph.get_endpoints()
    startpoints = graph.get_startpoints()
    translation_starts = graph.get_start_nodes()
    translation_stops = graph.get_stop_nodes()

    G.add_edges_from(edges.keys())
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
    
    edge_colors = []
    for edge in G.edges:
        if graph.edges[edges[edge]].edge_type == "translated":
            from_node = graph.edges[edges[edge]].from_node
            frame = graph.nodes[from_node].frame
        else:
            frame = None

        if frame == 0:
            edge_colors.append((1,0,0))
        elif frame == 1:
            edge_colors.append((0,1,0))
        elif frame == 2:
            edge_colors.append((0,0,1))
        else:
            edge_colors.append((0,0,0))
            


    nx.draw_networkx(G, pos=pos, node_shape='o', node_size=100, node_color=node_colors, edge_color=edge_colors, with_labels=False)

    plt.show()


if __name__ == "__main__":
    dg = RDG()
    dg.add_open_reading_frame(30, 90)
    dg.add_stop_codon_readthrough(90, 120)

    dg.add_open_reading_frame(131, 171)
    dg.add_open_reading_frame(150,850)

    plot(dg)
