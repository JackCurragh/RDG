import matplotlib.pyplot as plt
import networkx as nx
from RDG import RDG, Node, Edge
from matplotlib.gridspec import GridSpec



def position_graphs_nodes(graph) -> dict:
    '''
    determine the positioning of rdg nodes from graph structure. 
    This is a hacky method that determines Y axis by node key rather than calculating an optimum based on overlaps
    '''

    node_x_positions = [(graph.nodes[node].key, graph.nodes[node].node_start) for node in graph.nodes]
    pos = {node:(x_position, node) for node, x_position in node_x_positions}

    return pos


def position_horizontal_paths(graph) -> dict:
    '''
    Calculate node positions for RDG with horizontal paths. 
    '''
    pos = position_graphs_nodes(graph)
    for node in graph.nodes:
        # if node is a stop node, move it to the same y position (height) as the start node 
        if graph.nodes[node].node_type in ['stop', '3_prime']:
            upstream_node = graph.nodes[node].input_nodes[0]
            if graph.nodes[upstream_node].node_type in ['start', '5_prime', 'stop']:
                pos[node] = (pos[node][0], pos[upstream_node][1])
            else:
                pos[node] = (pos[node][0], pos[upstream_node][1] + 1) 

    # start the graph in line with the lowest numbered start node
    start_nodes = [node for node in graph.nodes if graph.nodes[node].node_type == 'start']
    pos[1] = (pos[1][0], min([pos[node][1] for node in start_nodes]))

    #hacky way to ensure that the non-coding path is at the top of the graph
    
    highest_node = max(list(pos.keys()))
    pos[2] = (pos[2][0], pos[highest_node][1] + 1)
    return pos

default_color_dict = {
    'edge_colors':{
        'frame0':(1,0,0),
        'frame1':(0,1,0),
        'frame2':(0,0,1)
    },
    'node_colors':{
        'startpoint':(0,0,1),
        'endpoint':(0.5,0,0.5),
        'translation_start':(0,1,0),
        'translation_stop':(1,0,0),
        'frameshift':(1,0.5,0.3)
    }
}




def plot(graph, color_dict=default_color_dict, node_size=50, edge_width=1.5, height_ratios=[2, 1], label_nodes=False, show_non_coding=False):
    # G = nx.DiGraph()
    G = nx.Graph()

    orfs_in_frame = {0:[], 1:[], 2:[]} 

    orfs = graph.get_orfs()
    for orf in orfs: 
        if len(orf) == 2:
            orfs_in_frame[orf[0]%3].append((orf[0], orf[1] - orf[0]))
        elif len(orf) == 3:
            orfs_in_frame[orf[0]%3].append((orf[0], orf[1] - orf[0]))
            orfs_in_frame[orf[1]%3].append((orf[1], orf[2] - orf[1]))



    pos = position_horizontal_paths(graph)

    endpoints = graph.get_endpoints()
    startpoints = graph.get_startpoints()
    translation_starts = graph.get_start_nodes()
    translation_stops = graph.get_stop_nodes()
    frameshifts = graph.get_frameshifts()
    print(translation_starts)

    if show_non_coding:
        pass
    else:
        graph.remove_edge(1)

    edges = graph.get_edges_from_to()

    #assign correct colouring based on frame to each ORF
    G.add_edges_from(edges.keys())
    node_colors = []
    for node in G.nodes():
        if node in startpoints:
            node_colors.append(color_dict['node_colors']['startpoint'])
        elif node in endpoints:
            node_colors.append(color_dict['node_colors']['endpoint'])
        elif node in translation_starts:
            node_colors.append(color_dict['node_colors']['translation_start'])
        elif node in translation_stops:
            node_colors.append(color_dict['node_colors']['translation_stop'])
        elif node in frameshifts:
            node_colors.append(color_dict['node_colors']['frameshift'])
        else:
            node_colors.append((0,0,0))
    
    edge_colors = []
    for edge in G.edges:
        if graph.edges[edges[edge]].edge_type == "translated":
            frame = graph.edges[edges[edge]].frame
        else:
            frame = None

        if frame == 0:
            edge_colors.append(color_dict['edge_colors']['frame0'])
        elif frame == 1:
            edge_colors.append(color_dict['edge_colors']['frame1'])
        elif frame == 2:
            edge_colors.append(color_dict['edge_colors']['frame2'])
        else:
            edge_colors.append((0,0,0))

    fig = plt.figure()
    fig.suptitle(f"Visualisation of an RDG for {graph.locus}")

    gs = GridSpec(2, 1, width_ratios=[1], height_ratios=height_ratios, hspace=0)
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])

    nx.draw_networkx(G, pos=pos, ax=ax1, node_shape='o', node_size=node_size, node_color=node_colors, width=edge_width, edge_color=edge_colors, with_labels=label_nodes)

    height = 10
    yticks_heights = []
    yticks_labels = []
    for frame, color in zip(sorted(orfs_in_frame.keys()), [(1,0,0), (0,1,0),  (0,0,1)]):
        yticks_heights.append(height + 1)
        yticks_labels.append("Frame " + str(frame + 1))
        ax2.broken_barh(sorted(orfs_in_frame[frame]), (height, 2), facecolors=color, edgecolor='black',)
        height += 3

    ax2.tick_params(bottom=True, labelbottom=True)
    ax2.set_xlim(left=ax1.get_xlim()[0], right=ax1.get_xlim()[1])
    ax2.set_ylim(bottom=10, top=19)
    ax2.set_yticks(yticks_heights, labels=yticks_labels)

    plt.show()
    return fig, ax1, ax2


if __name__ == "__main__":
    dg = RDG(name="Example Gene")
    dg.add_open_reading_frame(30, 90)
    dg.add_open_reading_frame(61, 400)
    dg.add_open_reading_frame(92, 150)
    dg.add_open_reading_frame(549, 849)

    no_node_color_dict = {
        'edge_colors':{
            'frame0':(1,0,0),
            'frame1':(0,1,0),
            'frame2':(0,0,1)
        },
        'node_colors':{
            'startpoint':(0,0,0),
            'endpoint':(0,0,0),
            'translation_start':(0,0,0),
            'translation_stop':(0,0,0),
            'frameshift':(0,0,0)
        }
    }


    plot(dg, color_dict=no_node_color_dict, edge_width=3, label_nodes=False)
