from RDG import Node

import unittest

def Node_types():
    node = Node(1, "error", (1,1))

class TestNode(unittest.TestCase):
        def test(self):
            with self.assertRaises(Exception) as context:
                Node_types()

            self.assertTrue('not valid. Valid types are' in str(context.exception))


