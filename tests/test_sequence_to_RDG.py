from RDG.sequence_to_RDG import extract_translons, build_graphs_from_fasta

import pytest


@pytest.fixture
def example_fasta_path():
    return "tests/test_data/test.fa"


def test_extract_translons():
    # Test case 1: Minimal input
    sequence = "ATGCTCTGA"
    result = extract_translons(sequence, min_length=1)
    assert result == [(0, 8)]

    # Test case 2: No start codon
    sequence = "CGTACCGTAGCTAG"
    result = extract_translons(sequence)
    assert result == []

    # Test case 3: Multiple start and stop codons
    sequence = "ATGTGTACCTAGTAGTATAA"
    result = extract_translons(sequence, starts=["ATG", "GTG"])
    assert result == [(0, 11), (2, 19)]

    # Test case 4: Custom start codons
    sequence = "CTGAAAAAAAAAAAATGCCCATAGTGA"
    result = extract_translons(sequence, starts=["CTG"])
    assert result == [(0, 23)]

    # Test case 5: Custom minimum length
    sequence = "ATGCGTACCTAGTAGTAGTAA"
    result = extract_translons(sequence, min_length=100)
    assert result == []


def test_build_graphs_from_fasta(example_fasta_path):
    # Test case 1: Minimal input
    result = build_graphs_from_fasta(
        example_fasta_path, start_codons=["ATG"], min_length=100
    )
    assert (
        result[0].newick()
        == "(((5:176)4:377,((8:64)7:332,((11:176)10:188,2:364)9:32)6:157)3:337)1;"
    )
