from RDG import RDG, Node, Edge
from sqlitedict import SqliteDict


def to_dict(obj):
    output = {}
    for key, item in obj.__dict__.items():
        if isinstance(item, list):
            l = []
            for item in item:
                if isinstance(item, int):
                    l.append(item)
                else:
                    d = to_dict(item)
                    l.append(d)
            output[key] = l
        else:
            output[key] = item

    return output


def save(graph, save_file):

    transcript = graph.locus
    graph_dict = {transcript: {}}
    for attr, value in graph.__dict__.items():
        if attr not in graph_dict[transcript]:
            graph_dict[transcript][attr] = {}

        if isinstance(value, dict):
            for key in value:
                graph_dict[transcript][attr][key] = to_dict(value[key])
        else:
            graph_dict[transcript][attr] = value

    try:
        with SqliteDict(save_file) as output_dict:
            for key in graph_dict:
                output_dict[key] = graph_dict[key]
                output_dict.commit()

    except Exception as ex:
        print("Error during storing data (Possibly unsupported):", ex)


def load(locus, cache_file="test_output.sqlite"):
    try:
        with SqliteDict(cache_file) as mydict:
            print(list(mydict.keys()))
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
