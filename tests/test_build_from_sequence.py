from RDG.build_RDG_from_sequence import parse_sequence_into_translated_regions, build_graphs_from_fasta

import pytest

@pytest.fixture
def example_fasta_path():
    return "tests/test_data/test.fa"

def test_parse_sequence_into_translated_regions():
    # Test case 1: Minimal input
    sequence = "ATGCTCTGA"
    result = parse_sequence_into_translated_regions(sequence, min_length=1)
    assert result == [(0, 9)]

    # Test case 2: No start codon
    sequence = "CGTACCGTAGCTAG"
    result = parse_sequence_into_translated_regions(sequence)
    assert result == []

    # Test case 3: Multiple start and stop codons
    sequence = "ATGTGTACCTAGTAGTATAA"
    result = parse_sequence_into_translated_regions(sequence, starts=["ATG", "GTG"])
    assert result == [(0, 12), (2, 20)]

    # Test case 4: Custom start codons
    sequence = "CTGAAAAAAAAAAAATGCCCATAGTGA"
    result = parse_sequence_into_translated_regions(sequence, starts=["CTG"])
    assert result == [(0, 24)]

    # Test case 5: Custom minimum length
    sequence = "ATGCGTACCTAGTAGTAGTAA"
    result = parse_sequence_into_translated_regions(sequence, min_length=100)
    assert result == []

def test_build_graphs_from_fasta(example_fasta_path):
    # Test case 1: Minimal input
    result = build_graphs_from_fasta(example_fasta_path, start_codons=["ATG"], min_lenth=100)
    assert result[0].newick() == '(((5:175)4:378,((8:63)7:333,((11:175)10:189,((14:63)13:201,2:264)12:100)9:32)6:157)3:337)1;'
