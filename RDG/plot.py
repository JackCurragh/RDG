import matplotlib.pyplot as plt
import matplotlib.patches as patches

from typing import Set

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
        branch_positions[branch] = graph.nodes[branch].position
        if branch not in upstream_branches:
            upstream_branches[branch] = [endpoint]
        else:
            upstream_branches[branch].append(endpoint)

    for branch in sorted(branch_positions, key=branch_positions.get,
                         reverse=True):
        for endpoint in upstream_branches[branch]:
            pos[endpoint] = (graph.nodes[endpoint].position, len(pos))

    while len(upstream_branches) >= 1:
        for branchpoint in upstream_branches:
            if len(upstream_branches[branchpoint]) == 2:
                # if the branchpoint has 2 downstream nodes then it is a
                # branchpoint and the nodes are placed at the same height
                pos[branchpoint] = (
                    graph.nodes[branchpoint].position,
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
            graph.nodes[startpoint].position, pos[out_node][1]
            )

    unassigned_nodes = [node for node in graph.nodes if node not in pos]
    for node in unassigned_nodes:
        upstream = graph.nodes[node].input_nodes[0]
        pos[node] = (graph.nodes[node].position, pos[upstream][1])

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


def get_reinitiation_nodes(graph, base_limit=0) -> (RDG, list):
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
            if stop_node.position > non_coding_edges[edge][0] + base_limit\
              and stop_node.position < non_coding_edges[edge][1]:
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


def score_kozak_context(sequence: str, position: int, window_size: int = 9) -> float:
    """
    Score the Kozak sequence context around a potential start codon.
    
    Kozak consensus: (gcc)gccRccAUGG (R = purine)
    Most important positions: -3 (should be purine) and +4 (should be G)
    
    Parameters:
    sequence (str): DNA sequence
    position (int): Position of start codon
    window_size (int): Size of context window to consider
    
    Returns:
    float: Kozak score (0-1, higher is better)
    """
    if position < window_size // 2 or position + window_size // 2 >= len(sequence):
        return 0.0, ""
    
    # Extract context around start codon
    context = sequence[position - window_size // 2:position + window_size // 2 + 1]
    
    if len(context) != window_size or not context[window_size // 2:window_size // 2 + 3] in ['ATG', 'CTG', 'GTG']:
        return 0.0, context
    
    score = 0.0
    
    # Position -3 (should be purine A or G)
    if position >= 3:
        if sequence[position - 3] in ['A', 'G']:
            score += 0.4
    
    # Position +4 (should be G)
    if position + 4 < len(sequence):
        if sequence[position + 4] == 'G':
            score += 0.4
    
    # Additional context scoring
    # Prefer G/C rich context upstream
    upstream = sequence[max(0, position - 6):position]
    gc_content = (upstream.count('G') + upstream.count('C')) / len(upstream) if upstream else 0
    score += 0.2 * gc_content
    
    return min(score, 1.0), context


def find_codons_with_context(sequence: str, 
                          start_codons: Set[str] = {"ATG", "CTG", "TTG", "ACG", "GTG"}):
    """
    Find ORFs with Kozak context scoring.
    
    Returns:
    List[Tuple[int, int, float]]: List of (start, stop, kozak_score) tuples
    """
    sequence = sequence.upper()
    stop_codons = {'TAA', 'TAG', 'TGA'}
    frames = {
        0: [],
        1: [],
        2: []
    }
    
    for frame in range(3):
        for i in range(frame, len(sequence) - 2, 3):
            codon = sequence[i:i+3]
            
            if codon in start_codons:
                kozak_score, context = score_kozak_context(sequence, i)
                frames[frame].append(
                    {
                        "codon": codon,
                        "kozak": kozak_score,
                        "context": context, 
                        "codon_type": "start",
                        "position": i
                    }
                )
               
            if codon in stop_codons:
                frames[frame].append(
                    {
                        "codon": codon,
                        "codon_type": "stop",
                        "position": i
                    }
                )     
    return frames

def plot_rdg(
    graph,
    ax,
    color_dict=default_color_dict,
    show_non_coding=True,
    translon_height=0.5,
    scantron_height=0.1,
    label_nodes=False,
    reinit_base_limit=0,
    allow_reinitiation=False,
):
    '''
    Plot the RDG decision graph on given axes

    :param graph: The RDG object representing the decision graph
    :type graph: RDG
    :param ax: The axes to plot on
    :type ax: matplotlib.axes.Axes
    :param color_dict: A dictionary of colors to use for the plot
    :type color_dict: dict
    :param show_non_coding: Whether to show the non coding edges
    :type show_non_coding: bool
    :param translon_height: The height of the translon plot
    :type translon_height: float
    :param scantron_height: The height of the non coding edges
    :type scantron_height: float
    :param label_nodes: Whether to label the nodes
    :type label_nodes: bool
    :param reinit_base_limit: The base limit for reinitiation
    :type reinit_base_limit: int
    :param allow_reinitiation: Whether to show reinitiation events
    :type allow_reinitiation: bool

    :return: The prepared graph object
    :rtype: RDG
    '''
    translons = []
    for translon in graph.get_translons():
        if translon not in translons:
            translons.append(translon)
    name = graph.name
    locus_stop = graph.locus_stop
    graph = RDG(name=name, locus_stop=locus_stop)
    for translon in translons:
        graph.add_open_reading_frame(translon[0], translon[1])

    # position nodes on xy plane
    pos = cladogram_layout(graph)

    if show_non_coding:
        pass
    else:
        graph.remove_edge(1)

    # Turn off axis labels
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel('')
    ax.set_ylabel('')

    max_x = max([pos[node][0] for node in pos])
    ax.set_xlim(0, max_x)

    max_y = max([pos[node][1] for node in pos])
    ax.set_ylim(0, max_y+2)

    # handle scaling of translon Y axis offset.
    translon_scaling = (translon_height / 2) - (scantron_height / 2)

    branch_heights = get_branch_heights(graph, pos)
    vertical_branch_width = graph.locus_stop * 0.01
    
    # Vertical lines at branch points
    for branch in graph.get_branch_points():
        start_node_coord = pos[branch]
        stop_node_coord = pos[graph.nodes[branch].output_nodes[0]]
        translon_length = stop_node_coord[0] - start_node_coord[0]

        vertical_branch_width = vertical_branch_width\
            if translon_length > vertical_branch_width else translon_length
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
            ax.add_patch(rect)

    # Only process reinitiation nodes if allow_reinitiation is True
    if allow_reinitiation:
        reinitiation_nodes = get_reinitiation_nodes(
            graph,
            base_limit=reinit_base_limit
            )
        # Vertical lines at reinitiation nodes
        for node in reinitiation_nodes:
            stop_node_coord = pos[node]
            from_node = graph.edges[reinitiation_nodes[node]].to_node

            reinitiation_edge_coord = (stop_node_coord[0], pos[from_node][1])
            base_height = min(stop_node_coord[1], reinitiation_edge_coord[1])

            out_node_y = pos[graph.nodes[node].output_nodes[0]][1]
            stop_node_height = out_node_y - translon_scaling

            height = abs(stop_node_height - reinitiation_edge_coord[1])

            width = vertical_branch_width/2

            start_node_coord = pos[graph.nodes[node].input_nodes[0]]
            translon_length = stop_node_coord[0] - start_node_coord[0]
            width = width if translon_length > width else translon_length
            rect = patches.Rectangle(
                (pos[node][0] - width, base_height),
                width=width,
                height=height,
                linewidth=0.5,
                edgecolor='#6d6d6d',
                facecolor='#6d6d6d',
                )
            ax.add_patch(rect)

    edges = graph.get_edges_from_to()
    # Draw edges
    for edge in edges:
        if edge[0] in pos and edge[1] in pos:
            if graph.edges[edges[edge]].edge_type == "translated":
                length = pos[edge[1]][0] - pos[edge[0]][0]

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
                ax.add_patch(rect)

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
                ax.add_patch(rect)

    if label_nodes:
        # label nodes in pos
        for node in pos:
            ax.text(
                pos[node][0],
                pos[node][1],
                node,
                horizontalalignment="center",
                verticalalignment="center",
                fontsize=18,
                color=color_dict["node_colors"]["startpoint"],
            )
    
    return graph


def plot_translon_frames(graph, ax, color_dict=default_color_dict):
    '''
    Plot translon frames on given axes
    
    :param graph: The RDG object representing the decision graph
    :type graph: RDG
    :param ax: The axes to plot on
    :type ax: matplotlib.axes.Axes
    :param color_dict: A dictionary of colors to use for the plot
    :type color_dict: dict
    
    :return: None
    '''
    # store translons in each frame for plotting
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
                (translon[1], translon[2] - translon[1])
                )
    
    # Plot translon frames
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
        ax.broken_barh(
            sorted(translons_in_frame[frame]),
            (height, 2),
            facecolors=color,
            edgecolor="black",
        )
        height += 3

    ax.tick_params(bottom=True, labelbottom=True)
    ax.set_ylim(bottom=10, top=19)
    ax.set_yticks(yticks_heights, labels=yticks_labels)
    ax.set_xlabel('Coordinates (nt)')


def plot_orf_structure(sequence, ax, color_dict=default_color_dict, max_coord=None):
    '''
    Plot ORF structure from sequence with codon annotations
    
    :param sequence: DNA sequence to analyze
    :type sequence: str
    :param ax: The axes to plot on
    :type ax: matplotlib.axes.Axes
    :param color_dict: A dictionary of colors to use for the plot
    :type color_dict: dict
    :param max_coord: Maximum coordinate for plot alignment
    :type max_coord: int
    
    :return: None
    '''
    frames = find_codons_with_context(sequence)
    
    height = 10
    yticks_heights = []
    yticks_labels = []
    
    for frame_num in [1, 2, 3]:
        frame_idx = frame_num - 1
        yticks_heights.append(height + 1)
        yticks_labels.append(f"Frame {frame_num}")
        
        frame_codons = frames.get(frame_idx, [])
        
        # Find ORFs in this frame (start to stop pairs)
        start_codons = [c for c in frame_codons if c["codon_type"] == "start"]
        stop_codons = [c for c in frame_codons if c["codon_type"] == "stop"]
        
        # Build list of significant start codons with their weights
        significant_starts = []
        for start_info in start_codons:
            start_pos = start_info["position"]
            start_codon = start_info.get("codon", "ATG")
            kozak_score = start_info.get("kozak", 0)
            
            if start_codon == "ATG" or kozak_score > 0.4:
                # Calculate ribosome weight: ATGs get base weight, then Kozak multiplies
                base_weight = 1.0 if start_codon == "ATG" else 0.5  # ATG stronger than near-cognates
                ribosome_weight = base_weight * (0.3 + 0.7 * kozak_score)  # Kozak modulates strength
                significant_starts.append((start_pos, ribosome_weight, start_codon))
        
        # Calculate ribosome density at each position
        sequence_length = int(max_coord) if max_coord else len(sequence)
        ribosome_density = [0.0] * sequence_length
        
        for start_pos, weight, start_codon in significant_starts:
            # Find next stop codon in same frame
            next_stop = None
            for stop_info in stop_codons:
                stop_pos = stop_info["position"]
                if stop_pos > start_pos:
                    next_stop = stop_pos
                    break
            
            # Add ribosome density from this start to its stop
            if next_stop:
                for pos in range(start_pos, min(next_stop, sequence_length)):
                    ribosome_density[pos] += weight
        
        # Plot ORFs with density-based coloring
        plotted_regions = set()  # Track plotted regions to avoid overlap
        
        for start_pos, weight, start_codon in significant_starts:
            # Find stop for this start
            next_stop = None
            for stop_info in stop_codons:
                stop_pos = stop_info["position"]
                if stop_pos > start_pos:
                    next_stop = stop_pos
                    break
            
            if next_stop:
                # Plot ORF region with density-based intensity
                region_key = (start_pos, next_stop)
                if region_key not in plotted_regions:
                    plotted_regions.add(region_key)
                    
                    # Get average density for this region
                    avg_density = sum(ribosome_density[start_pos:next_stop]) / (next_stop - start_pos)
                    
                    # Map density to green intensity (0-3+ density -> light to dark green)
                    if avg_density > 2.0:
                        orf_color = "#004d00"  # Very dark green for high ribosome load
                    elif avg_density > 1.5:
                        orf_color = "#006600"  # Dark green
                    elif avg_density > 1.0:
                        orf_color = "#008000"  # Medium-dark green
                    elif avg_density > 0.5:
                        orf_color = "#32cd32"  # Medium green
                    else:
                        orf_color = "#90ee90"  # Light green for low density
                    
                    # Plot ORF as continuous block
                    orf_rect = patches.Rectangle(
                        (start_pos, height),
                        next_stop - start_pos,  # Full ORF length
                        2,  # full frame height
                        facecolor=orf_color,
                        alpha=1.0,  # Solid color to show density clearly
                        edgecolor="black",
                        linewidth=0.5,
                        zorder=5
                    )
                    ax.add_patch(orf_rect)
        
        # Mark individual start codons with markers
        for start_pos, weight, start_codon in significant_starts:
            bar_color = "#000080" if start_codon == "ATG" else "#800080"  # Blue for ATG, purple for others
            start_rect = patches.Rectangle(
                (start_pos, height),
                8,  # Thick start marker
                2,  # full frame height
                facecolor=bar_color,
                alpha=1.0,
                edgecolor=bar_color,  # Same color as bar for clean look
                linewidth=0,
                zorder=10
            )
            ax.add_patch(start_rect)
            
        # Plot frame background - distinct colors to clearly show frames
        frame_colors = ["#f0f0f0", "#e0e0e0", "#d0d0d0"]  # Light gray gradients  
        frame_color = frame_colors[frame_idx]
        background_rect = patches.Rectangle(
            (0, height),
            max_coord or len(sequence),
            2,
            facecolor=frame_color,
            alpha=1.0,  # Solid background for clear frame distinction
            edgecolor="black",
            linewidth=1,
            zorder=1  # Put background behind everything
        )
        ax.add_patch(background_rect)
        
        height += 3
    
    ax.tick_params(bottom=True, labelbottom=True)
    ax.set_ylim(bottom=10, top=19)
    ax.set_yticks(yticks_heights, labels=yticks_labels)
    ax.set_xlabel('Coordinates (nt)')
    
    if max_coord:
        ax.set_xlim(0, max_coord)
    else:
        ax.set_xlim(0, len(sequence))
    


def compose(
    graph,
    plot_type="translons",
    sequence=None,
    color_dict=default_color_dict,
    height_ratios=[2, 1],
    show_non_coding=True,
    translon_height=0.5,
    scantron_height=0.1,
    label_nodes=False,
    reinit_base_limit=0,
    allow_reinitiation=False,
):
    '''
    Compose RDG visualization with different bottom panel options
    
    :param graph: The RDG object representing the decision graph
    :type graph: RDG
    :param plot_type: Type of bottom plot ("translons" or "orfs")
    :type plot_type: str
    :param sequence: DNA sequence for ORF analysis (required if plot_type="orfs")
    :type sequence: str
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
    :param label_nodes: Whether to label the nodes
    :type label_nodes: bool
    :param reinit_base_limit: The base limit for reinitiation
    :type reinit_base_limit: int
    :param allow_reinitiation: Whether to show reinitiation events
    :type allow_reinitiation: bool
    
    :return: The figure and axes objects
    :rtype: tuple
    '''
    # Set up the figure
    fig = plt.figure()
    fig.suptitle(f"Visualisation of an RDG for {graph.name}")

    gs = GridSpec(2, 1, width_ratios=[1],
                  height_ratios=height_ratios, hspace=0)
    ax1 = fig.add_subplot(gs[0])  # RDG plot
    ax2 = fig.add_subplot(gs[1])  # Bottom plot

    # Plot the RDG
    processed_graph = plot_rdg(
        graph,
        ax1,
        color_dict=color_dict,
        show_non_coding=show_non_coding,
        translon_height=translon_height,
        scantron_height=scantron_height,
        label_nodes=label_nodes,
        reinit_base_limit=reinit_base_limit,
        allow_reinitiation=allow_reinitiation,
    )
    
    # Align x-axes
    ax2.set_xlim(left=ax1.get_xlim()[0], right=ax1.get_xlim()[1])
    
    # Plot bottom panel based on type
    if plot_type == "translons":
        plot_translon_frames(processed_graph, ax2, color_dict=color_dict)
    elif plot_type == "orfs":
        if sequence is None:
            raise ValueError("Sequence required for ORF structure plot")
        plot_orf_structure(
            sequence, 
            ax2, 
            color_dict=color_dict, 
            max_coord=ax1.get_xlim()[1]
        )
    else:
        raise ValueError(f"Unknown plot_type: {plot_type}. Use 'translons' or 'orfs'")
    
    return fig, ax1, ax2


# For backward compatibility
def plot(
    graph,
    color_dict=default_color_dict,
    height_ratios=[2, 1],
    show_non_coding=True,
    translon_height=0.5,
    scantron_height=0.1,
    label_nodes=False,
    reinit_base_limit=0,
    allow_reinitiation=False,
    sequence=None,
):
    '''
    Legacy plot function - calls compose with translons plot type
    '''
    return compose(
        graph,
        plot_type="translons",
        sequence=sequence,
        color_dict=color_dict,
        height_ratios=height_ratios,
        show_non_coding=show_non_coding,
        translon_height=translon_height,
        scantron_height=scantron_height,
        label_nodes=label_nodes,
        reinit_base_limit=reinit_base_limit,
        allow_reinitiation=allow_reinitiation,
    )


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
    # g.add_open_readin node_color_dict)
    g = RDG(name="ATF4 - NM_001675", locus_stop=2041)
    g.add_open_reading_frame(486, 649)
    g.add_open_reading_frame(700, 891, reinitiation=False)
    g.add_open_reading_frame(888, 1943, reinitiation=False)

    plot(g, color_dict=no_node_color_dict, reinit_base_limit=300)
    plot(g, color_dict=no_node_color_dict, reinit_base_limit=0)

    # g = RDG(name="test", locus_stop=2041)
    # g.add_open_reading_frame(486, 649)
    # g.add_open_reading_frame(700, 891)
    # g.add_open_reading_frame(888, 1943)
    # g.add_open_reading_frame(85, 94)
    # g.add_open_reading_frame(23, 26)

    # plot(g, color_dict=no_node_color_dict)
