from RDG import RDG, Node, Edge, save, load

# from RDG_to_file import save, load
import unittest


def test_read_write():
    dg = RDG()
    dg.add_open_reading_frame(30, 90)
    dg.add_open_reading_frame(131, 171)
    dg.add_open_reading_frame(150, 850)

    save(dg, "test_output.sqlite")
    dg2 = load("name", "test_output.sqlite")

    assert dg.describe() == dg2.describe()


def invalid_load_file():
    dg = RDG()
    dg.add_open_reading_frame(30, 90)
    dg.add_open_reading_frame(131, 171)
    dg.add_open_reading_frame(150, 850)

    save(dg, "test_output.sqlite")
    dg2 = load("name", "")


class TestAddSCR(unittest.TestCase):
    def test(self):
        with self.assertRaises(Exception) as context:
            invalid_load_file()

        self.assertTrue("Error during loading data" in str(context.exception))
