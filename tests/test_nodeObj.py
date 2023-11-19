from RDG import RDG, Node

import unittest


def Node_types_error():
    node = Node(1, "error", (1, 1))


class TestNode(unittest.TestCase):
    def test(self):
        with self.assertRaises(Exception) as context:
            Node_types_error()

        self.assertTrue("not valid. Valid types are" in str(context.exception))


def test_node_init():
    g = RDG()
    node_key = g.get_new_node_key()
    node = Node(
        key=node_key,
        node_type="stop",
        position=5,
        edges_in=[],
        edges_out=[],
        nodes_in=[],
        nodes_out=[],
    )
    assert node.key == node_key
    assert node.node_type == "stop"
    assert node.node_start == 5


def test_node_key():
    g = RDG()
    node_key = g.get_new_node_key()
    node = Node(
        key=node_key,
        node_type="stop",
        position=5,
        edges_in=[],
        edges_out=[],
        nodes_in=[],
        nodes_out=[],
    )
    assert node.key == node_key


def test_node_type():
    g = RDG()
    node_key = g.get_new_node_key()
    node = Node(
        key=node_key,
        node_type="stop",
        position=5,
        edges_in=[],
        edges_out=[],
        nodes_in=[],
        nodes_out=[],
    )
    assert node.node_type == "stop"
