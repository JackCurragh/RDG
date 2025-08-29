from RDG import RDG, Node

import pytest


def Node_types_error():
    node = Node(1, "error", (1, 1))


def test_node_types_error():
    with pytest.raises(ValueError, match="Invalid node type: error"):
        Node_types_error()


def test_node_init():
    g = RDG()
    node_key = g.get_new_node_key()
    node = Node(
        key=node_key,
        node_type="stop",
        position=5,
        input_edges=[],
        output_edges=[],
        input_nodes=[],
        output_nodes=[],
    )
    assert node.key == node_key
    assert node.node_type == "stop"
    assert node.position == 5


def test_node_key():
    g = RDG()
    node_key = g.get_new_node_key()
    node = Node(
        key=node_key,
        node_type="stop",
        position=5,
        input_edges=[],
        output_edges=[],
        input_nodes=[],
        output_nodes=[],
    )
    assert node.key == node_key


def test_node_type():
    g = RDG()
    node_key = g.get_new_node_key()
    node = Node(
        key=node_key,
        node_type="stop",
        position=5,
        input_edges=[],
        output_edges=[],
        input_nodes=[],
        output_nodes=[],
    )
    assert node.node_type == "stop"
