import matplotlib.pyplot as plt
import networkx as nx
from RDG import RDG, Node, Edge
from matplotlib.gridspec import GridSpec



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
    frameshifts = graph.get_frameshifts()

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
        elif node in frameshifts:
            node_colors.append((1,0.5,0.3))
        else:
            node_colors.append((0,0,0))
    
    edge_colors = []
    for edge in G.edges:
        if graph.edges[edges[edge]].edge_type == "translated":
            frame = graph.edges[edges[edge]].frame
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
            
    fig = plt.figure()
    fig.suptitle("Visualisation of an RDG")

    gs = GridSpec(2, 1, width_ratios=[1], height_ratios=[4, 1], hspace=0)
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])

    nx.draw_networkx(G, pos=pos, ax=ax1, node_shape='o', node_size=100, node_color=node_colors, edge_color=edge_colors, with_labels=False)

    orfs_in_frame = {0:[], 1:[], 2:[]} 

    orfs = graph.get_orfs()
    for orf in orfs: 
        if len(orf) == 2:
            orfs_in_frame[orf[0]%3].append((orf[0], orf[1] - orf[0]))
        elif len(orf) == 3:
            orfs_in_frame[orf[0]%3].append((orf[0], orf[1] - orf[0]))
            orfs_in_frame[orf[1]%3].append((orf[1], orf[2] - orf[1]))



    height = 10
    yticks_heights = []
    yticks_labels = []
    for frame, color in zip(sorted(orfs_in_frame.keys()), [(1,0,0), (0,1,0),  (0,0,1)]):
        yticks_heights.append(height + 1)
        yticks_labels.append("Frame " + str(frame + 1))
        ax2.broken_barh(sorted(orfs_in_frame[frame])[::-1], (height, 2), facecolors=color, edgecolor='black',)
        height += 3

    ax2.tick_params(bottom=True, labelbottom=True)
    ax2.set_xlim(left=ax1.get_xlim()[0], right=ax1.get_xlim()[1])
    ax2.set_ylim(bottom=10, top=19)
    ax2.set_yticks(yticks_heights, labels=yticks_labels)

    plt.show()


if __name__ == "__main__":
    dg = RDG()
    dg.add_open_reading_frame(30, 90)
    dg.add_open_reading_frame(61, 400)
    dg.add_open_reading_frame(92, 150)

    # dg.add_open_reading_frame(550, 850)
    # dg.add_stop_codon_readthrough(850, 880)
    dg.add_frameshift(400, 450, 2)
    # print(dg.get_orfs())
    plot(dg)
