import consensus
from consensus import AlleleParser, ComparisonResult, ConsensusAlgorithm


def test_comparison_results():
    a1 = AlleleParser("prog1", "A*02:01")
    a2 = AlleleParser("prog1", "A*02:01:01:01")
    assert a1.is_equal(a2) == ComparisonResult.MORE_RESOLUTION
    assert a2.is_equal(a1) == ComparisonResult.LESS_RESOLUTION
    assert a1.is_equal(a1) == ComparisonResult.EQUAL

    a3 = AlleleParser("prog1", "B*02:01:01:01")
    assert a1.is_equal(a3) == ComparisonResult.NOT_EQUAL


def test_allele_parsing():
    a1 = AlleleParser("prog1", "A*02:01")
    str(a1)
    a1 = AlleleParser("prog1", "A*02")
    str(a1)
    a1 = AlleleParser("prog1", "")
    str(a1)


def test_consensus():
    a1 = AlleleParser("prog1", "A*02:01")
    a2 = AlleleParser("prog1", "A*02:01:01:01")
    con = ConsensusAlgorithm(a1.gene)
    con.add(a1)
    con.add(a2)
    assert con.alleles[a2] == 2


def test_overlapping_allele_consensus():
    allele_list = [
        "DPA1*01 (0.6666)",
        "DPA1*02 (0.3335)",
        "DPA1*01:03:01 (0.3333)",
        "DPA1*01:04 (0.3333)",
        "DPA1*01:03:01:01 (0.3114)",
    ]
    parsed_list = [AlleleParser("prog1", a) for a in allele_list]
    con = ConsensusAlgorithm("DPA1")
    for a in parsed_list:
        con.add(a)
    assert con.alleles[parsed_list[3]] == 2
    assert con.alleles[parsed_list[4]] == 3
