from RDG import RDG, Node, Edge
from sqlitedict import SqliteDict


def to_dict(obj):
    '''
    Recursively convert a class to a dictionary

    Parameters
    ----------
    obj : class
        The class to convert to a dictionary

    Returns 
    -------
    output : dict
    '''
    output = {}
    for key, item in obj.__dict__.items():
        if isinstance(item, list):
            l = []
            for item in item:
                l.append(item)

            output[key] = l
        else:
            output[key] = item

    return output


def save(graph, save_file):
    '''
    Save a graph to a file

    Parameters
    ----------
    graph : RDG

    save_file : str
        The name of the file to save to (should end in .sqlite)
    '''
    graph_dict = {graph.locus: {}}
    for attr, value in graph.__dict__.items():
        if attr not in graph_dict[graph.locus]:
            graph_dict[graph.locus][attr] = {}

        if isinstance(value, dict):
            for key in value:
                graph_dict[graph.locus][attr][key] = to_dict(value[key])
        else:
            graph_dict[graph.locus][attr] = value

    try:
        with SqliteDict(save_file) as output_dict:
            for key in graph_dict:
                output_dict[key] = graph_dict[key]
                output_dict.commit()

    except:
        raise Exception("Error during storing data (Possibly unsupported):")
    
def newick_to_file(newick, save_file):
    '''
    Save a newick string to a file

    Parameters
    ----------
    newick : str
        The newick string to save

    save_file : str
        The name of the file to save to (should end in .sqlite)
    '''
    with open(save_file, 'w') as file:
        file.write(newick)





def load(locus, cache_file="test_output.sqlite") -> RDG:
    '''
    Load a graph from a file

    Parameters
    ----------
    locus : str
        The locus of the graph to load

    cache_file : str
        The name of the file to load from (should end in .sqlite)
    '''
    try:
        with SqliteDict(cache_file) as mydict:
            graph_dict = mydict[
                locus
            ]  # No need to use commit(), since we are only loading data!
    except Exception as ex:
        raise Exception("Error during loading data:", ex)

    dg = RDG()
    nodes = {}
    for key in graph_dict["nodes"]:
        nodes[key] = Node(
            key,
            graph_dict["nodes"][key]["node_type"],
            graph_dict["nodes"][key]["node_start"],
            graph_dict["nodes"][key]["input_edges"],
            graph_dict["nodes"][key]["output_edges"],
            graph_dict["nodes"][key]["input_nodes"],
            graph_dict["nodes"][key]["output_nodes"],
        )

    edges = {}

    for key in graph_dict["edges"]:
        edges[key] = Edge(
            key,
            graph_dict["edges"][key]["edge_type"],
            graph_dict["edges"][key]["from_node"],
            graph_dict["edges"][key]["to_node"],
            graph_dict["edges"][key]["coordinates"],
        )

    dg = dg.Load(
        locus_name=locus,
        locus_start=graph_dict["locus_start"],
        locus_stop=graph_dict["locus_stop"],
        nodes=nodes,
        edges=edges,
    )
    return dg
