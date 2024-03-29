import pytest

from index import Gene


@pytest.fixture()
def original_sequence():
    return (
        "AGGGGAAGTGACATAGCAGGAACTACTAGTACCCTTCAGGAACAAATAGGATGGATGACACATAATCCACCTATCC"
        "CAGTAGGAGAAATCTATAAAAGATGGATAATCCTGGGATTAAATAAAATAGTAAGAATGTATAGCCCTACCAGCAT"
        "TCTGGACATAAGACAAGGACCAAAGGAACCCTTTAGAGACTATGTAGACCGATTCTATAAAACTCTAAGAGCCGAG"
        "CAAGCTTCACAAGAGGTAAAAAATTGGATGACAGAAACCTTGTTGGTCCAAAATGCGAACCCAGATTGTAAGACTA"
        "TTTTAAAAGCATTGGGACCAGGAGCGACACTAGAAGAAATGATGACAGCATGTCAGGGAGTGGGGGGACCCGGCCA"
        "TAAAGCAAGAGTTTTGGCTGAAGCAATGAGCCAAGTAACAAATCCAGCTACCATAATGATACAGAAAGGCAATTTT"
        "AGGAACCAAAGAAAGACTGTTAAGTGTTTCAATTGTGGCAAAGAAGGGCACATAGCC"
    )


def test_minimum_gap(original_sequence):
    gaps_to_test = [0, 12, 20, 100]

    for gap in gaps_to_test:
        gene = Gene(
            original_sequence=original_sequence, minimum_CpG_gap=gap, gap_method=1
        )
        if gap == 0:
            assert gene.new_sequence == (
                "CGAGGATCGGACATCGCGGGAACGACGTCGACGCTTCAGGAACAAATCGGATGGATGACG"
                "CATAATCCGCCGATCCCGGTCGGCGAAATCTATAAACGATGGATAATCCTCGGATTAAAT"
                "AAAATCGTACGAATGTATTCGCCGACGTCGATTCTCGACATACGACAAGGACCGAAGGAA"
                "CCGTTTCGAGACTACGTCGACCGATTCTATAAAACGCTACGAGCCGAGCAAGCGTCGCAA"
                "GAGGTAAAAAATTGGATGACGGAAACGTTGCTCGTCCAAAACGCGAACCCGGATTGTAAG"
                "ACGATTTTAAAAGCGCTCGGACCGGGCGCGACGCTCGAAGAAATGATGACGGCGTGTCAG"
                "GGCGTCGGCGGACCCGGCCATAAAGCGCGAGTTCTCGCGGAAGCGATGTCGCAAGTAACG"
                "AATCCGGCGACGATAATGATACAGAAAGGCAATTTTCGAAACCAACGAAAGACGGTTAAG"
                "TGTTTCAATTGCGGCAAAGAAGGGCACATCGCG"
            )
        elif gap == 12:
            assert gene.new_sequence == (
                "CGAGGAAGTGACATCGCAGGAACTACTTCGACCCTTCAGGAACAAATCGGATGGATGACA"
                "CATAATCCGCCTATCCCAGTAGGCGAAATCTATAAAAGATGGATAATCCTCGGATTAAAT"
                "AAAATCGTAAGAATGTATTCGCCTACCAGCATTCTCGACATAAGACAAGGACCGAAGGAA"
                "CCCTTTCGAGACTATGTAGACCGATTCTATAAAACTCTAAGAGCCGAGCAAGCTTCACAA"
                "GAGGTAAAAAATTGGATGACGGAAACCTTGTTGGTCCAAAATGCGAACCCAGATTGTAAG"
                "ACGATTTTAAAAGCATTGGGACCAGGAGCGACACTAGAAGAAATGATGACGGCATGTCAG"
                "GGAGTGGGGGGACCCGGCCATAAAGCAAGAGTTCTCGCTGAAGCAATGTCGCAAGTAACA"
                "AATCCGGCTACCATAATGATACAGAAAGGCAATTTTCGAAACCAAAGAAAGACGGTTAAG"
                "TGTTTCAATTGCGGCAAAGAAGGGCACATCGCC"
            )
        elif gap == 20:
            assert gene.new_sequence == (
                "CGAGGAAGTGACATAGCAGGAACGACTAGTACCCTTCAGGAACAAATCGGATGGATGACA"
                "CATAATCCACCGATCCCAGTAGGAGAAATCTATAAACGATGGATAATCCTGGGATTAAAT"
                "AAAATCGTAAGAATGTATAGCCCTACCTCGATTCTGGACATAAGACAAGGACCGAAGGAA"
                "CCCTTTAGAGACTATGTAGACCGATTCTATAAAACTCTAAGAGCCGAGCAAGCTTCACAA"
                "GAGGTAAAAAATTGGATGACGGAAACCTTGTTGGTCCAAAATGCGAACCCAGATTGTAAG"
                "ACTATTTTAAAAGCATTGGGACCAGGAGCGACACTAGAAGAAATGATGACAGCGTGTCAG"
                "GGAGTGGGGGGACCCGGCCATAAAGCAAGAGTTTTGGCGGAAGCAATGAGCCAAGTAACA"
                "AATCCGGCTACCATAATGATACAGAAAGGCAATTTTCGAAACCAAAGAAAGACTGTTAAG"
                "TGTTTCAATTGCGGCAAAGAAGGGCACATAGCC"
            )
        elif gap == 100:
            assert gene.new_sequence == (
                "CGAGGAAGTGACATAGCAGGAACTACTAGTACCCTTCAGGAACAAATAGGATGGATGACA"
                "CATAATCCACCTATCCCAGTAGGAGAAATCTATAAAAGATGGATAATCCTGGGATTAAAT"
                "AAAATAGTAAGAATGTATAGCCCTACCAGCATTCTGGACATAAGACAAGGACCAAAGGAA"
                "CCCTTTAGAGACTATGTAGACCGATTCTATAAAACTCTAAGAGCCGAGCAAGCTTCACAA"
                "GAGGTAAAAAATTGGATGACAGAAACCTTGTTGGTCCAAAATGCGAACCCAGATTGTAAG"
                "ACTATTTTAAAAGCATTGGGACCAGGAGCGACACTAGAAGAAATGATGACAGCATGTCAG"
                "GGAGTGGGGGGACCCGGCCATAAAGCAAGAGTTTTGGCTGAAGCAATGAGCCAAGTAACA"
                "AATCCAGCTACCATAATGATACAGAAAGGCAATTTTAGGAACCAAAGAAAGACTGTTAAG"
                "TGTTTCAATTGCGGCAAAGAAGGGCACATAGCC"
            )


# Testing as above but with added A enrichment
def test_A_rich(original_sequence):
    gaps_to_test = [0, 12, 20, 100]

    for gap in gaps_to_test:
        gene = Gene(
            original_sequence=original_sequence, minimum_CpG_gap=gap, gap_method=1
        )
        gene.mutate_A_rich()
        if gap == 0:
            assert gene.new_sequence == (
                "CGAGGATCGGACATCGCGGGAACGACGTCGACGCTACAAGAACAAATCGGATGGATGACG"
                "CATAATCCGCCGATACCGGTCGGCGAAATATATAAACGATGGATAATACTCGGATTAAAT"
                "AAAATCGTACGAATGTATTCGCCGACGTCGATACTCGACATACGACAAGGACCGAAGGAA"
                "CCGTTTCGAGACTACGTCGACCGATTCTATAAAACGCTACGAGCCGAACAAGCGTCGCAA"
                "GAAGTAAAAAATTGGATGACGGAAACGTTACTCGTACAAAACGCGAACCCGGATTGTAAG"
                "ACGATATTAAAAGCGCTCGGACCGGGCGCGACGCTCGAAGAAATGATGACGGCGTGTCAA"
                "GGCGTCGGCGGACCCGGACATAAAGCGCGAGTACTCGCGGAAGCGATGTCGCAAGTAACG"
                "AATCCGGCGACGATAATGATACAAAAAGGAAATTTTCGAAACCAACGAAAGACGGTAAAG"
                "TGTTTCAATTGCGGAAAAGAAGGACACATCGCG"
            )
        elif gap == 12:
            assert gene.new_sequence == (
                "CGAGGAAGTGACATCGCAGGAACAACATCGACACTACAAGAACAAATCGGATGGATGACA"
                "CATAATCCGCCAATACCAGTAGGCGAAATATATAAAAGATGGATAATACTCGGATTAAAT"
                "AAAATCGTAAGAATGTATTCGCCAACAAGCATACTCGACATAAGACAAGGACCGAAGGAA"
                "CCATTTCGAGACTATGTAGACCGATTCTATAAAACACTAAGAGCCGAACAAGCATCACAA"
                "GAAGTAAAAAATTGGATGACGGAAACATTATTAGTACAAAATGCGAACCCAGATTGTAAG"
                "ACGATATTAAAAGCATTAGGACCAGGAGCGACACTAGAAGAAATGATGACGGCATGTCAA"
                "GGAGTAGGAGGACCCGGACATAAAGCAAGAGTACTCGCAGAAGCAATGTCGCAAGTAACA"
                "AATCCGGCAACAATAATGATACAAAAAGGAAATTTTCGAAACCAAAGAAAGACGGTAAAG"
                "TGTTTCAATTGCGGAAAAGAAGGACACATCGCA"
            )
        elif gap == 20:
            assert gene.new_sequence == (
                "CGAGGAAGTGACATAGCAGGAACGACAAGTACACTACAAGAACAAATCGGATGGATGACA"
                "CATAATCCACCGATACCAGTAGGAGAAATATATAAACGATGGATAATACTAGGATTAAAT"
                "AAAATCGTAAGAATGTATAGCCCAACATCGATACTAGACATAAGACAAGGACCGAAGGAA"
                "CCATTTAGAGACTATGTAGACCGATTCTATAAAACACTAAGAGCCGAACAAGCATCACAA"
                "GAAGTAAAAAATTGGATGACGGAAACATTATTAGTACAAAATGCGAACCCAGATTGTAAG"
                "ACAATATTAAAAGCATTAGGACCAGGAGCGACACTAGAAGAAATGATGACAGCGTGTCAA"
                "GGAGTAGGAGGACCCGGACATAAAGCAAGAGTATTAGCGGAAGCAATGAGCCAAGTAACA"
                "AATCCGGCAACAATAATGATACAAAAAGGAAATTTTCGAAACCAAAGAAAGACAGTAAAG"
                "TGTTTCAATTGCGGAAAAGAAGGACACATAGCA"
            )
        elif gap == 100:
            assert gene.new_sequence == (
                "CGAGGAAGTGACATAGCAGGAACAACAAGTACACTACAAGAACAAATAGGATGGATGACA"
                "CATAATCCACCAATACCAGTAGGAGAAATATATAAAAGATGGATAATACTAGGATTAAAT"
                "AAAATAGTAAGAATGTATAGCCCAACAAGCATACTAGACATAAGACAAGGACCAAAGGAA"
                "CCATTTAGAGACTATGTAGACCGATTCTATAAAACACTAAGAGCCGAACAAGCATCACAA"
                "GAAGTAAAAAATTGGATGACAGAAACATTATTAGTACAAAATGCGAACCCAGATTGTAAG"
                "ACAATATTAAAAGCATTAGGACCAGGAGCGACACTAGAAGAAATGATGACAGCATGTCAA"
                "GGAGTAGGAGGACCCGGACATAAAGCAAGAGTATTAGCAGAAGCAATGAGCCAAGTAACA"
                "AATCCAGCAACAATAATGATACAAAAAGGAAATTTTAGAAACCAAAGAAAGACAGTAAAG"
                "TGTTTCAATTGCGGAAAAGAAGGACACATAGCA"
            )


# Testing with minimum gap of 0 as it creates the most CpG's
def test_packaging_signal_protection(original_sequence):
    # 55 tests a previous bug
    packaging_signal_lengths_to_test = [1, 2, 3, 4, 55, 100]

    for length in packaging_signal_lengths_to_test:
        gene = Gene(
            original_sequence=original_sequence,
            minimum_CpG_gap=0,
            gap_method=1,
            packaging_signal_length_beginning=length,
            packaging_signal_length_end=length,
        )
        if length == 1:
            assert gene.new_sequence == (
                "AGGGGATCGGACATCGCGGGAACGACGTCGACGCTTCAGGAACAAATCGGATGGATGACG"
                "CATAATCCGCCGATCCCGGTCGGCGAAATCTATAAACGATGGATAATCCTCGGATTAAAT"
                "AAAATCGTACGAATGTATTCGCCGACGTCGATTCTCGACATACGACAAGGACCGAAGGAA"
                "CCGTTTCGAGACTACGTCGACCGATTCTATAAAACGCTACGAGCCGAGCAAGCGTCGCAA"
                "GAGGTAAAAAATTGGATGACGGAAACGTTGCTCGTCCAAAACGCGAACCCGGATTGTAAG"
                "ACGATTTTAAAAGCGCTCGGACCGGGCGCGACGCTCGAAGAAATGATGACGGCGTGTCAG"
                "GGCGTCGGCGGACCCGGCCATAAAGCGCGAGTTCTCGCGGAAGCGATGTCGCAAGTAACG"
                "AATCCGGCGACGATAATGATACAGAAAGGCAATTTTCGAAACCAACGAAAGACGGTTAAG"
                "TGTTTCAATTGCGGCAAAGAAGGGCACATCGCC"
            )
        elif length == 2:
            assert gene.new_sequence == (
                "AGGGGATCGGACATCGCGGGAACGACGTCGACGCTTCAGGAACAAATCGGATGGATGACG"
                "CATAATCCGCCGATCCCGGTCGGCGAAATCTATAAACGATGGATAATCCTCGGATTAAAT"
                "AAAATCGTACGAATGTATTCGCCGACGTCGATTCTCGACATACGACAAGGACCGAAGGAA"
                "CCGTTTCGAGACTACGTCGACCGATTCTATAAAACGCTACGAGCCGAGCAAGCGTCGCAA"
                "GAGGTAAAAAATTGGATGACGGAAACGTTGCTCGTCCAAAACGCGAACCCGGATTGTAAG"
                "ACGATTTTAAAAGCGCTCGGACCGGGCGCGACGCTCGAAGAAATGATGACGGCGTGTCAG"
                "GGCGTCGGCGGACCCGGCCATAAAGCGCGAGTTCTCGCGGAAGCGATGTCGCAAGTAACG"
                "AATCCGGCGACGATAATGATACAGAAAGGCAATTTTCGAAACCAACGAAAGACGGTTAAG"
                "TGTTTCAATTGCGGCAAAGAAGGGCACATCGCC"
            )
        elif length == 3:
            assert gene.new_sequence == (
                "AGGGGATCGGACATCGCGGGAACGACGTCGACGCTTCAGGAACAAATCGGATGGATGACG"
                "CATAATCCGCCGATCCCGGTCGGCGAAATCTATAAACGATGGATAATCCTCGGATTAAAT"
                "AAAATCGTACGAATGTATTCGCCGACGTCGATTCTCGACATACGACAAGGACCGAAGGAA"
                "CCGTTTCGAGACTACGTCGACCGATTCTATAAAACGCTACGAGCCGAGCAAGCGTCGCAA"
                "GAGGTAAAAAATTGGATGACGGAAACGTTGCTCGTCCAAAACGCGAACCCGGATTGTAAG"
                "ACGATTTTAAAAGCGCTCGGACCGGGCGCGACGCTCGAAGAAATGATGACGGCGTGTCAG"
                "GGCGTCGGCGGACCCGGCCATAAAGCGCGAGTTCTCGCGGAAGCGATGTCGCAAGTAACG"
                "AATCCGGCGACGATAATGATACAGAAAGGCAATTTTCGAAACCAACGAAAGACGGTTAAG"
                "TGTTTCAATTGCGGCAAAGAAGGGCACATAGCC"
            )
        elif length == 4:
            assert gene.new_sequence == (
                "AGGGGATCGGACATCGCGGGAACGACGTCGACGCTTCAGGAACAAATCGGATGGATGACG"
                "CATAATCCGCCGATCCCGGTCGGCGAAATCTATAAACGATGGATAATCCTCGGATTAAAT"
                "AAAATCGTACGAATGTATTCGCCGACGTCGATTCTCGACATACGACAAGGACCGAAGGAA"
                "CCGTTTCGAGACTACGTCGACCGATTCTATAAAACGCTACGAGCCGAGCAAGCGTCGCAA"
                "GAGGTAAAAAATTGGATGACGGAAACGTTGCTCGTCCAAAACGCGAACCCGGATTGTAAG"
                "ACGATTTTAAAAGCGCTCGGACCGGGCGCGACGCTCGAAGAAATGATGACGGCGTGTCAG"
                "GGCGTCGGCGGACCCGGCCATAAAGCGCGAGTTCTCGCGGAAGCGATGTCGCAAGTAACG"
                "AATCCGGCGACGATAATGATACAGAAAGGCAATTTTCGAAACCAACGAAAGACGGTTAAG"
                "TGTTTCAATTGCGGCAAAGAAGGGCACATAGCC"
            )
        elif length == 55:
            assert gene.new_sequence == (
                "AGGGGAAGTGACATAGCAGGAACTACTAGTACCCTTCAGGAACAAATAGGATGGATGACG"
                "CATAATCCGCCGATCCCGGTCGGCGAAATCTATAAACGATGGATAATCCTCGGATTAAAT"
                "AAAATCGTACGAATGTATTCGCCGACGTCGATTCTCGACATACGACAAGGACCGAAGGAA"
                "CCGTTTCGAGACTACGTCGACCGATTCTATAAAACGCTACGAGCCGAGCAAGCGTCGCAA"
                "GAGGTAAAAAATTGGATGACGGAAACGTTGCTCGTCCAAAACGCGAACCCGGATTGTAAG"
                "ACGATTTTAAAAGCGCTCGGACCGGGCGCGACGCTCGAAGAAATGATGACGGCGTGTCAG"
                "GGCGTCGGCGGACCCGGCCATAAAGCGCGAGTTCTCGCGGAAGCGATGTCGCAAGTAACG"
                "AATCCGGCGACGATAATGATACAGAAAGGCAATTTTAGGAACCAAAGAAAGACTGTTAAG"
                "TGTTTCAATTGTGGCAAAGAAGGGCACATAGCC"
            )
        elif length == 100:
            assert gene.new_sequence == (
                "AGGGGAAGTGACATAGCAGGAACTACTAGTACCCTTCAGGAACAAATAGGATGGATGACA"
                "CATAATCCACCTATCCCAGTAGGAGAAATCTATAAAAGATGGATAATCCTCGGATTAAAT"
                "AAAATCGTACGAATGTATTCGCCGACGTCGATTCTCGACATACGACAAGGACCGAAGGAA"
                "CCGTTTCGAGACTACGTCGACCGATTCTATAAAACGCTACGAGCCGAGCAAGCGTCGCAA"
                "GAGGTAAAAAATTGGATGACGGAAACGTTGCTCGTCCAAAACGCGAACCCGGATTGTAAG"
                "ACGATTTTAAAAGCGCTCGGACCGGGCGCGACGCTCGAAGAAATGATGACGGCGTGTCAG"
                "GGCGTCGGCGGACCCGGCCATAAAGCGCGAGTTCTCGCGGAAGCGATGTCGCAAGTAACA"
                "AATCCAGCTACCATAATGATACAGAAAGGCAATTTTAGGAACCAAAGAAAGACTGTTAAG"
                "TGTTTCAATTGTGGCAAAGAAGGGCACATAGCC"
            )

    # Also check when making A-rich
    for length in packaging_signal_lengths_to_test:
        gene = Gene(
            original_sequence=original_sequence,
            minimum_CpG_gap=0,
            gap_method=1,
            packaging_signal_length_beginning=length,
            packaging_signal_length_end=length,
        )
        gene.mutate_A_rich()
        if length == 1:
            assert gene.new_sequence == (
                "AGAGGATCGGACATCGCGGGAACGACGTCGACGCTACAAGAACAAATCGGATGGATGACG"
                "CATAATCCGCCGATACCGGTCGGCGAAATATATAAACGATGGATAATACTCGGATTAAAT"
                "AAAATCGTACGAATGTATTCGCCGACGTCGATACTCGACATACGACAAGGACCGAAGGAA"
                "CCGTTTCGAGACTACGTCGACCGATTCTATAAAACGCTACGAGCCGAACAAGCGTCGCAA"
                "GAAGTAAAAAATTGGATGACGGAAACGTTACTCGTACAAAACGCGAACCCGGATTGTAAG"
                "ACGATATTAAAAGCGCTCGGACCGGGCGCGACGCTCGAAGAAATGATGACGGCGTGTCAA"
                "GGCGTCGGCGGACCCGGACATAAAGCGCGAGTACTCGCGGAAGCGATGTCGCAAGTAACG"
                "AATCCGGCGACGATAATGATACAAAAAGGAAATTTTCGAAACCAACGAAAGACGGTAAAG"
                "TGTTTCAATTGCGGAAAAGAAGGACACATCGCC"
            )
        elif length == 2:
            assert gene.new_sequence == (
                "AGAGGATCGGACATCGCGGGAACGACGTCGACGCTACAAGAACAAATCGGATGGATGACG"
                "CATAATCCGCCGATACCGGTCGGCGAAATATATAAACGATGGATAATACTCGGATTAAAT"
                "AAAATCGTACGAATGTATTCGCCGACGTCGATACTCGACATACGACAAGGACCGAAGGAA"
                "CCGTTTCGAGACTACGTCGACCGATTCTATAAAACGCTACGAGCCGAACAAGCGTCGCAA"
                "GAAGTAAAAAATTGGATGACGGAAACGTTACTCGTACAAAACGCGAACCCGGATTGTAAG"
                "ACGATATTAAAAGCGCTCGGACCGGGCGCGACGCTCGAAGAAATGATGACGGCGTGTCAA"
                "GGCGTCGGCGGACCCGGACATAAAGCGCGAGTACTCGCGGAAGCGATGTCGCAAGTAACG"
                "AATCCGGCGACGATAATGATACAAAAAGGAAATTTTCGAAACCAACGAAAGACGGTAAAG"
                "TGTTTCAATTGCGGAAAAGAAGGACACATCGCC"
            )
        elif length == 3:
            assert gene.new_sequence == (
                "AGGGGATCGGACATCGCGGGAACGACGTCGACGCTACAAGAACAAATCGGATGGATGACG"
                "CATAATCCGCCGATACCGGTCGGCGAAATATATAAACGATGGATAATACTCGGATTAAAT"
                "AAAATCGTACGAATGTATTCGCCGACGTCGATACTCGACATACGACAAGGACCGAAGGAA"
                "CCGTTTCGAGACTACGTCGACCGATTCTATAAAACGCTACGAGCCGAACAAGCGTCGCAA"
                "GAAGTAAAAAATTGGATGACGGAAACGTTACTCGTACAAAACGCGAACCCGGATTGTAAG"
                "ACGATATTAAAAGCGCTCGGACCGGGCGCGACGCTCGAAGAAATGATGACGGCGTGTCAA"
                "GGCGTCGGCGGACCCGGACATAAAGCGCGAGTACTCGCGGAAGCGATGTCGCAAGTAACG"
                "AATCCGGCGACGATAATGATACAAAAAGGAAATTTTCGAAACCAACGAAAGACGGTAAAG"
                "TGTTTCAATTGCGGAAAAGAAGGACACATAGCC"
            )
        elif length == 4:
            assert gene.new_sequence == (
                "AGGGGATCGGACATCGCGGGAACGACGTCGACGCTACAAGAACAAATCGGATGGATGACG"
                "CATAATCCGCCGATACCGGTCGGCGAAATATATAAACGATGGATAATACTCGGATTAAAT"
                "AAAATCGTACGAATGTATTCGCCGACGTCGATACTCGACATACGACAAGGACCGAAGGAA"
                "CCGTTTCGAGACTACGTCGACCGATTCTATAAAACGCTACGAGCCGAACAAGCGTCGCAA"
                "GAAGTAAAAAATTGGATGACGGAAACGTTACTCGTACAAAACGCGAACCCGGATTGTAAG"
                "ACGATATTAAAAGCGCTCGGACCGGGCGCGACGCTCGAAGAAATGATGACGGCGTGTCAA"
                "GGCGTCGGCGGACCCGGACATAAAGCGCGAGTACTCGCGGAAGCGATGTCGCAAGTAACG"
                "AATCCGGCGACGATAATGATACAAAAAGGAAATTTTCGAAACCAACGAAAGACGGTAAAG"
                "TGTTTCAATTGCGGAAAAGAAGGACACATAGCC"
            )
        elif length == 55:
            assert gene.new_sequence == (
                "AGGGGAAGTGACATAGCAGGAACTACTAGTACCCTTCAGGAACAAATAGGATGGATGACG"
                "CATAATCCGCCGATACCGGTCGGCGAAATATATAAACGATGGATAATACTCGGATTAAAT"
                "AAAATCGTACGAATGTATTCGCCGACGTCGATACTCGACATACGACAAGGACCGAAGGAA"
                "CCGTTTCGAGACTACGTCGACCGATTCTATAAAACGCTACGAGCCGAACAAGCGTCGCAA"
                "GAAGTAAAAAATTGGATGACGGAAACGTTACTCGTACAAAACGCGAACCCGGATTGTAAG"
                "ACGATATTAAAAGCGCTCGGACCGGGCGCGACGCTCGAAGAAATGATGACGGCGTGTCAA"
                "GGCGTCGGCGGACCCGGACATAAAGCGCGAGTACTCGCGGAAGCGATGTCGCAAGTAACG"
                "AATCCGGCGACGATAATGATACAAAAAGGAAATTTTAGGAACCAAAGAAAGACTGTTAAG"
                "TGTTTCAATTGTGGCAAAGAAGGGCACATAGCC"
            )
        elif length == 100:
            assert gene.new_sequence == (
                "AGGGGAAGTGACATAGCAGGAACTACTAGTACCCTTCAGGAACAAATAGGATGGATGACA"
                "CATAATCCACCTATCCCAGTAGGAGAAATCTATAAAAGATGGATAATACTCGGATTAAAT"
                "AAAATCGTACGAATGTATTCGCCGACGTCGATACTCGACATACGACAAGGACCGAAGGAA"
                "CCGTTTCGAGACTACGTCGACCGATTCTATAAAACGCTACGAGCCGAACAAGCGTCGCAA"
                "GAAGTAAAAAATTGGATGACGGAAACGTTACTCGTACAAAACGCGAACCCGGATTGTAAG"
                "ACGATATTAAAAGCGCTCGGACCGGGCGCGACGCTCGAAGAAATGATGACGGCGTGTCAA"
                "GGCGTCGGCGGACCCGGACATAAAGCGCGAGTACTCGCGGAAGCGATGTCGCAAGTAACA"
                "AATCCAGCTACCATAATGATACAGAAAGGCAATTTTAGGAACCAAAGAAAGACTGTTAAG"
                "TGTTTCAATTGTGGCAAAGAAGGGCACATAGCC"
            )


def test_optimizing_average_gap(original_sequence):
    gene = Gene(
        original_sequence=original_sequence,
        gap_method=2,
        minimum_CpG_gap=None,
        desired_CpG_gap=29,
    )
    assert gene.minimum_CpG_gap == 20

    # Same as in test_minimum_gap
    assert gene.new_sequence == (
        "CGAGGAAGTGACATAGCAGGAACGACTAGTACCCTTCAGGAACAAATCGGATGGATGACA"
        "CATAATCCACCGATCCCAGTAGGAGAAATCTATAAACGATGGATAATCCTGGGATTAAAT"
        "AAAATCGTAAGAATGTATAGCCCTACCTCGATTCTGGACATAAGACAAGGACCGAAGGAA"
        "CCCTTTAGAGACTATGTAGACCGATTCTATAAAACTCTAAGAGCCGAGCAAGCTTCACAA"
        "GAGGTAAAAAATTGGATGACGGAAACCTTGTTGGTCCAAAATGCGAACCCAGATTGTAAG"
        "ACTATTTTAAAAGCATTGGGACCAGGAGCGACACTAGAAGAAATGATGACAGCGTGTCAG"
        "GGAGTGGGGGGACCCGGCCATAAAGCAAGAGTTTTGGCGGAAGCAATGAGCCAAGTAACA"
        "AATCCGGCTACCATAATGATACAGAAAGGCAATTTTCGAAACCAAAGAAAGACTGTTAAG"
        "TGTTTCAATTGCGGCAAAGAAGGGCACATAGCC"
    )
