import matplotlib.pyplot as plt
import networkx as nx
from RDG import RDG, Node, Edge
from RDG_to_file import save, load
from matplotlib.gridspec import GridSpec

from ete3 import Tree


def layout_graph(graph: RDG, branch_height=500) -> dict:
    """
    Method for laying out the graph based on branch points 

    """
    pos = {}
    
    paths = graph.get_unique_paths()
    node_y_positions = {paths[0][0]: 0}

    branch_points = graph.get_branch_points()

    branch_num = 0 
    for path in paths:
        for i, node in enumerate(path):
            if node in branch_points:
                if node not in node_y_positions:
                    node_y_positions[node] = branch_num
                    branch_num += branch_height
            
            else:
                if node not in node_y_positions:
                    if graph.nodes[node].node_type == "3_prime":
                        upstream_node = graph.nodes[node].input_nodes[0]
                        if upstream_node in branch_points:
                            node_y_positions[node] = node_y_positions[path[i-1]] + 3
                        else:
                            node_y_positions[node] = node_y_positions[path[i-1]]
                    else:
                        node_y_positions[node] = node_y_positions[path[i-1]]

    for node in graph.nodes:
        pos[node] = (graph.nodes[node].node_start, node_y_positions[node])

    return pos



default_color_dict = {
    "edge_colors": {"frame0": (1, 0, 0), "frame1": (0, 1, 0), "frame2": (0, 0, 1)},
    "node_colors": {
        "startpoint": (0, 0, 1),
        "endpoint": (0.5, 0, 0.5),
        "translation_start": (0, 1, 0),
        "translation_stop": (1, 0, 0),
        "frameshift": (1, 0.5, 0.3),
    },
}


def get_tree_positioning(newick):
    """
    plot graph data strucuture using ete3 package
    """

    t = Tree(newick, format=1)
    t.render("mytree.png", w=183, units="mm")


def plot(
    graph,
    color_dict=default_color_dict,
    node_size=50,
    edge_width=1.5,
    height_ratios=[2, 1],
    label_nodes=False,
    show_non_coding=False,
):
    
    # The networkx graph must either be directed or the edges must be sorted before adding to graph
    G = nx.DiGraph() 



    # store orfs in each frame for plotting the ORF plot.
    orfs_in_frame = {0: [], 1: [], 2: []}
    orfs = graph.get_orfs()
    for orf in orfs:
        if len(orf) == 2:
            orfs_in_frame[orf[0] % 3].append((orf[0], orf[1] - orf[0]))
        elif len(orf) == 3:
            orfs_in_frame[orf[0] % 3].append((orf[0], orf[1] - orf[0]))
            orfs_in_frame[orf[1] % 3].append((orf[1], orf[2] - orf[1]))

    # position nodes on xy plane
    pos = layout_graph(graph)


    # identify feature types for colouring
    endpoints = graph.get_endpoints()
    startpoints = graph.get_startpoints()
    translation_starts = graph.get_start_nodes()
    translation_stops = graph.get_stop_nodes()
    frameshifts = graph.get_frameshifts()

    if show_non_coding:
        pass
    else:
        graph.remove_edge(1)
    
    
    edges = graph.get_edges_from_to()

    # assign correct colouring based on frame to each ORF
    G.add_edges_from(edges.keys())

    node_colors = []
    for node in G.nodes():
        if node in startpoints:
            node_colors.append(color_dict["node_colors"]["startpoint"])
        elif node in endpoints:
            node_colors.append(color_dict["node_colors"]["endpoint"])
        elif node in translation_starts:
            node_colors.append(color_dict["node_colors"]["translation_start"])
        elif node in translation_stops:
            node_colors.append(color_dict["node_colors"]["translation_stop"])
        elif node in frameshifts:
            node_colors.append(color_dict["node_colors"]["frameshift"])
        else:
            node_colors.append((0, 0, 0))

    edge_colors = []

    for edge in G.edges:
        if graph.edges[edges[edge]].edge_type == "translated":
            frame = graph.edges[edges[edge]].get_frame()
        else:
            frame = None

        if frame == 0:
            edge_colors.append(color_dict["edge_colors"]["frame0"])
        elif frame == 1:
            edge_colors.append(color_dict["edge_colors"]["frame1"])
        elif frame == 2:
            edge_colors.append(color_dict["edge_colors"]["frame2"])
        else:
            edge_colors.append((0, 0, 0))

    # set up the figure
    fig = plt.figure()
    fig.suptitle(f"Visualisation of an RDG for {graph.locus}")

    gs = GridSpec(2, 1, width_ratios=[1], height_ratios=height_ratios, hspace=0)
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])

    # plot the graph
    nx.draw_networkx(
        G,
        pos=pos,
        ax=ax1,
        node_shape="o",
        node_size=node_size,
        node_color=node_colors,
        width=edge_width,
        edge_color=edge_colors,
        with_labels=label_nodes,
        alpha=0.7,
        font_size=15,
    )

    # plot the orf plot
    height = 10
    yticks_heights = []
    yticks_labels = []
    for frame, color in zip(
        sorted(orfs_in_frame.keys()), [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
    ):
        yticks_heights.append(height + 1)
        yticks_labels.append("Frame " + str(frame + 1))
        ax2.broken_barh(
            sorted(orfs_in_frame[frame]),
            (height, 2),
            facecolors=color,
            edgecolor="black",
        )
        height += 3

    ax2.tick_params(bottom=True, labelbottom=True)
    ax2.set_xlim(left=ax1.get_xlim()[0], right=ax1.get_xlim()[1])
    ax2.set_ylim(bottom=10, top=19)
    ax2.set_yticks(yticks_heights, labels=yticks_labels)

    plt.show()
    return fig, ax1, ax2


if __name__ == "__main__":
    no_node_color_dict = {
        "edge_colors": {"frame0": (1, 0, 0), "frame1": (0, 1, 0), "frame2": (0, 0, 1)},
        "node_colors": {
            "startpoint": (0, 0, 0),
            "endpoint": (0, 0, 0),
            "translation_start": (0, 0, 0),
            "translation_stop": (0, 0, 0),
            "frameshift": (0, 0, 0),
        },
    }

    dg = RDG()
    dg.add_open_reading_frame(100, 200, reinitiation=True, upstream_limit=1)
    dg.add_open_reading_frame(300, 500, reinitiation=True, upstream_limit=1)
    dg.add_open_reading_frame(600, 700, reinitiation=True, upstream_limit=1)
    dg.add_open_reading_frame(800, 900, reinitiation=True, upstream_limit=1)

    print()
    plot(dg, color_dict=no_node_color_dict, edge_width=3, label_nodes=False, show_non_coding=True)

    # # print(dg.newick())
    # get_tree_positioning(dg.newick())

    # plot(dg, color_dict=no_node_color_dict, edge_width=3, label_nodes=True, show_non_coding=True)
