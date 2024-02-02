import matplotlib.pyplot as plt
import matplotlib.patches as patches

from RDG import RDG

# from RDG_to_file import save, load, newick_to_file
from matplotlib.gridspec import GridSpec


def cladogram_layout(graph: RDG) -> dict:
    """
    Method for the cladogram inspired RDG layout

    :param graph: The RDG object representing the decision graph
    :type graph: RDG

    :return: A dictionary representing the layout of the cladogram
    :rtype: dict
    """
    pos = {}  # node: (x, y)

    # endpoints define the height of the graph
    # The order of the endpoints ius determined by the upstream branch points
    # of the endpoints
    # The first endpoint is the one with the furthest 3' upstream branch point
    endpoints = graph.get_endpoints()
    upstream_branches = {}
    branch_positions = {}
    for endpoint in endpoints:
        branch = graph.get_upstream_branchpoint(endpoint)
        branch_positions[branch] = graph.nodes[branch].node_start
        if branch not in upstream_branches:
            upstream_branches[branch] = [endpoint]
        else:
            upstream_branches[branch].append(endpoint)

    for branch in sorted(branch_positions, key=branch_positions.get,
                         reverse=True):
        for endpoint in upstream_branches[branch]:
            pos[endpoint] = (graph.nodes[endpoint].node_start, len(pos))

    while len(upstream_branches) >= 1:
        for branchpoint in upstream_branches:
            if len(upstream_branches[branchpoint]) == 2:
                # if the branchpoint has 2 downstream nodes then it is a
                # branchpoint and the nodes are placed at the same height
                pos[branchpoint] = (
                    graph.nodes[branchpoint].node_start,
                    (
                        pos[upstream_branches[branchpoint][0]][1]
                        + pos[upstream_branches[branchpoint][1]][1]
                        ) / 2,
                )
                upstream_branch = graph.get_upstream_branchpoint(branchpoint)
                if upstream_branch not in upstream_branches:
                    upstream_branches[upstream_branch] = [branchpoint]
                else:
                    upstream_branches[upstream_branch].append(branchpoint)
                del upstream_branches[branchpoint]
                break

        if branchpoint in graph.get_startpoints():
            break

    startpoints = graph.get_startpoints()
    for startpoint in startpoints:
        out_node = graph.nodes[startpoint].output_nodes[0]
        pos[startpoint] = (
            graph.nodes[startpoint].node_start, pos[out_node][1]
            )

    unassigned_nodes = [node for node in graph.nodes if node not in pos]
    for node in unassigned_nodes:
        upstream = graph.nodes[node].input_nodes[0]
        pos[node] = (graph.nodes[node].node_start, pos[upstream][1])

    return pos


# used in the calculation of branch heights this function looks downstream of
# a node and returns the first end or branch node it finds.
# if the node itself is a branch or end node it returns itself as this is run
# on the downstream nodes of a branch point already
def get_end_or_branch(graph, node):
    '''
    Get the end or branch node immediately downstream of a node

    :param graph: The RDG object representing the decision graph
    :type graph: RDG
    :param node: The node to get the end or branch node of
    :type node: int

    :return: The end or branch node of the node
    :rtype: int
    '''
    if node in graph.get_endpoints() or node in graph.get_branch_points():
        return node
    for downstream_node in graph.nodes[node].output_nodes:
        if downstream_node in graph.get_endpoints() \
          or downstream_node in graph.get_branch_points():
            return downstream_node
        else:
            return get_end_or_branch(graph, downstream_node)


# Once node positions are know this function can be used to
# calculate the heights of the vertical branches that can be
# used to connect translons and non coding edges
def get_branch_heights(graph, pos):
    '''
    calculate the heights of the vertical branches that connect
    translons and non coding edges at branch points

    :param graph: The RDG object representing the decision graph
    :type graph: RDG

    :return: A dictionary of branch heights
    :rtype: dict
    '''
    branch_heights = {}
    for branch in graph.get_branch_points():
        A = get_end_or_branch(graph, graph.nodes[branch].output_nodes[0])
        B = get_end_or_branch(graph, graph.nodes[branch].output_nodes[1])
        branch_heights[branch] = abs(pos[A][1] - pos[B][1])
    return branch_heights


def get_reinitiation_nodes(graph) -> (RDG, list):
    '''
    identify reinitiation nodes and return the edges that they should
    reinitiate translation at

    :param graph: The RDG object representing the decision graph
    :type graph: RDG

    :return: A dict of reinitiation nodes
    :rtype: dict
    '''
    reinitiation_nodes = {}
    non_coding_edges = {
        edge: graph.edges[edge].coordinates for edge in graph.edges
        if graph.edges[edge].edge_type == "untranslated" and
        graph.edges[edge].to_node not in graph.get_endpoints()
            }

    for node in graph.get_stop_nodes():
        stop_node = graph.nodes[node]
        for edge in non_coding_edges:
            if stop_node.node_start > non_coding_edges[edge][0]\
              and stop_node.node_start < non_coding_edges[edge][1]:
                if edge != 1:  # ignore the non-coding path
                    reinitiation_nodes[node] = edge
                break

    return reinitiation_nodes


default_color_dict = {
        "edge_colors": {
            0: "#ffbb8d",
            1: "#ffeedd",
            2: "#ffd8be"
            },
        "node_colors": {
            "startpoint": "#000000",
            "endpoint": "#000000",
            "translation_start": "#00b050",
            "translation_stop": "#000000",
            "frameshift": "#000000",
        },
    }


def plot(
    graph,
    color_dict=default_color_dict,
    height_ratios=[2, 1],
    show_non_coding=True,
    translon_height=0.5,
    scantron_height=0.1,
    label_nodes=False,
):
    '''
    Generate a plot of the RDG

    :param graph: The RDG object representing the decision graph
    :type graph: RDGot

    :param color_dict: A dictionary of colors to use for the plot
    :type color_dict: dict

    :param height_ratios: The height ratios of the two subplots
    :type height_ratios: list

    :param show_non_coding: Whether to show the non coding edges
    :type show_non_coding: bool

    :param translon_height: The height of the translon plot
    :type translon_height: float

    :param scantron_height: The height of the non coding edges
    :type scantron_height: float

    :return: The figure and axes objects
    :rtype: tuple
    '''
    translons = []
    for translon in graph.get_translons():
        if translon not in translons:
            translons.append(translon)
    name = graph.locus
    locus_stop = graph.locus_stop
    graph = RDG(name=name, locus_stop=locus_stop)
    for translon in translons:
        graph.add_open_reading_frame(translon[0], translon[1])

    reinitiation_nodes = get_reinitiation_nodes(graph)
    # store translons in each frame for plotting the translon plot.
    translons_in_frame = {0: [], 1: [], 2: []}
    translons = graph.get_translons()
    for translon in translons:
        if len(translon) == 2:
            translons_in_frame[translon[0] % 3].append(
                (translon[0], translon[1] - translon[0])
                )
        elif len(translon) == 3:
            translons_in_frame[translon[0] % 3].append(
                (translon[0], translon[1] - translon[0])
                )
            translons_in_frame[translon[1] % 3].append(
                (translon[1], translon[2] - translon[1])Systematic analysis of the PTEN 5’ leader identifies a major AUU initiated proteoform
                )

    # position nodes on xy plane
    pos = cladogram_layout(graph)

    if show_non_coding:
        pass
    else:
        graph.remove_edge(1)

    # set up the figure
    fig = plt.figure()
    fig.suptitle(f"Visualisation of an RDG for {graph.locus}")

    gs = GridSpec(2, 1, width_ratios=[1],
                  height_ratios=height_ratios, hspace=0)
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])

    # Turn off axis labels
    ax1.set_xticks([])
    ax1.set_yticks([])
    ax1.set_xlabel('')
    ax1.set_ylabel('')

    max_x = max([pos[node][0] for node in pos])
    ax1.set_xlim(0, max_x)
    ax2.set_xlim(0, max_x)

    max_y = max([pos[node][1] for node in pos])
    ax1.set_ylim(0, max_y+2)
    ax2.set_ylim(0, max_y)

    # handle scaling of translon Y axis offset.
    # If translons are plotted without offset they do not align
    # with the non coding edges
    translon_scaling = (translon_height / 2) - (scantron_height / 2)

    branch_heights = get_branch_heights(graph, pos)
    vertical_branch_width = graph.locus_stop * 0.01
    # Vertical lines at branch points
    for branch in graph.get_branch_points():
        if branch in pos:
            A = get_end_or_branch(graph, graph.nodes[branch].output_nodes[0])
            B = get_end_or_branch(graph, graph.nodes[branch].output_nodes[1])
            base_height = min(pos[A][1], pos[B][1])
            rect = patches.Rectangle(
                (pos[branch][0], base_height),
                width=vertical_branch_width,
                height=branch_heights[branch],
                linewidth=0.5,
                edgecolor='#000000',
                facecolor='#000000',
                )
            ax1.add_patch(rect)

    # Vertical lines at reinitiation nodes
    for node in reinitiation_nodes:
        stop_node_coord = pos[node]
        from_node = graph.edges[reinitiation_nodes[node]].to_node

        reinitiation_edge_coord = (stop_node_coord[0], pos[from_node][1])
        base_height = min(stop_node_coord[1], reinitiation_edge_coord[1])

        out_node_y = pos[graph.nodes[node].output_nodes[0]][1]
        stop_node_height = out_node_y - translon_scaling

        height = abs(stop_node_height - reinitiation_edge_coord[1])

        width = min(
            vertical_branch_width/2,
            abs(reinitiation_edge_coord[0] - stop_node_coord[0])
        )
        rect = patches.Rectangle(
            (pos[node][0] - width, base_height),
            width=width,
            height=height,
            linewidth=0.5,
            edgecolor='#6d6d6d',
            facecolor='#6d6d6d',
            )
        ax1.add_patch(rect)

    edges = graph.get_edges_from_to()
    # loop over edges and make translated edges 10x thicker than
    # non-translated edges

    for edge in edges:
        if edge[0] in pos and edge[1] in pos:
            if graph.edges[edges[edge]].edge_type == "translated":
                length = pos[edge[1]][0] - pos[edge[0]][0]

                # Translon needs to be positioned at the same height as the
                # next node downstream and centered using the scaling
                ds_node_y = pos[graph.nodes[edge[1]].output_nodes[0]][1]
                frame = graph.nodes[edge[0]].frame

                rect = patches.Rectangle(
                    (pos[edge[0]][0], ds_node_y - translon_scaling),
                    length,
                    translon_height,
                    linewidth=0.5,
                    edgecolor='#000000',
                    facecolor=color_dict["edge_colors"][frame],
                    )
                ax1.add_patch(rect)

            else:
                length = pos[edge[1]][0] - pos[edge[0]][0]
                rect = patches.Rectangle(
                    (pos[edge[0]][0], pos[edge[1]][1]),
                    length,
                    scantron_height,
                    linewidth=0.5,
                    edgecolor='#000000',
                    facecolor='#000000',
                    )
                ax1.add_patch(rect)

    if label_nodes:
        # label nodes in pos
        for node in pos:
            ax1.text(
                pos[node][0],
                pos[node][1],
                node,
                horizontalalignment="center",
                verticalalignment="center",
                fontsize=18,
                color=color_dict["node_colors"]["startpoint"],
            )
    # plot the translon plot
    height = 10
    yticks_heights = []
    yticks_labels = []
    for frame, color in zip(
        sorted(translons_in_frame.keys()),
        [
            color_dict["edge_colors"][0],
            color_dict["edge_colors"][1],
            color_dict["edge_colors"][2]
            ]
    ):
        yticks_heights.append(height + 1)
        yticks_labels.append("Frame " + str(frame + 1))
        ax2.broken_barh(
            sorted(translons_in_frame[frame]),
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
        "edge_colors": {
            0: "#ffbb8d",
            1: "#ffeedd",
            2: "#ffd8be"
            },
        "node_colors": {
            "startpoint": "#000000",
            "endpoint": "#000000",
            "translation_start": "#00b050",
            "translation_stop": "#000000",
            "frameshift": "#000000",
        },
    }
    # g = RDG(name="SRD5A1 - NM_001047")
    # g.add_open_reading_frame(138, 917)
    # g.add_open_reading_frame(86, 517)
    # g.add_open_reading_frame(112, 438)
    # plot(g, color_dict=no_node_color_dict)

    # g = RDG(name="ATF4 - NM_001675", locus_stop=2041)
    # g.add_open_reading_frame(200, 293)
    # g.add_open_reading_frame(486, 1943)
    # g.add_open_reading_frame(700, 891, reinitiation=False)
    # g.add_open_reading_frame(888, 1943, reinitiation=False)

    # plot(g, color_dict=no_node_color_dict)

    g = RDG(name="test")
    g.add_open_reading_frame(10, 100)
    g.add_open_reading_frame(110, 200)
    g.add_open_reading_frame(210, 300)
    g.add_open_reading_frame(310, 400)
    g.add_open_reading_frame(410, 500)

    plot(g, color_dict=no_node_color_dict)
