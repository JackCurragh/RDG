from RDG import RDG, Node, Edge, save, load

# from RDG_to_file import save, load


def test_read_write():
    dg = RDG()
    dg.add_open_reading_frame(30, 90)
    dg.add_open_reading_frame(131, 171)
    dg.add_open_reading_frame(150, 850)

    save(dg, "test_output.sqlite")
    dg2 = load("", "test_output.sqlite")

    assert dg.describe() == dg2.describe()
