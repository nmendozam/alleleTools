from .consensus import (Allele, CmpResult, ConsensusAlgorithm, Report,
                        get_allele_pair)


def make_report(alleles: list):
    calls = {
        "sample": "",
        "calls": {
            "": {
                "": alleles,
            }
        },
    }
    return calls


def test_comparison_results():
    a1 = Allele("A*02:01")
    a2 = Allele("A*02:01:01:01")

    assert a1.compare(a2) == CmpResult.MORE_RESOLUTION
    assert a2.compare(a1) == CmpResult.LESS_RESOLUTION
    assert a1.compare(a1) == CmpResult.EQUAL

    a3 = Allele("B*02:01:01:01")
    assert a1.compare(a3) == CmpResult.NOT_EQUAL


def test_allele_parsing():
    a1 = Allele("A*02:01")
    assert str(a1) == "A*02:01"

    a1 = Allele("A*02")
    assert str(a1) == "A*02"


def test_allele_parsing_no_allele():
    failed = False
    try:
        a1 = Allele("")
    except:
        failed = True
    assert failed


def test_consensus():
    report = [
        {
            "gene": "A",
            "program": "alg1",
            "allele": Allele("A*02:01"),
        },
        {
            "gene": "A",
            "program": "alg2",
            "allele": Allele("A*02:01:01:01"),
        }
    ]
    con = ConsensusAlgorithm(report)
    assert len(con.get_flat_alleles()) == 2


def test_overlapping_allele_consensus():
    allele_list = [
        "DPA1*01 (0.6666)",
        "DPA1*02 (0.3335)",
        "DPA1*01:03:01 (0.3333)",
        "DPA1*01:04 (0.3333)",
        "DPA1*01:03:01:01 (0.3114)",
    ]
    calls = make_report(allele_list)
    report = Report(calls, resolution=2)
    con = ConsensusAlgorithm(report)
    print(con.get_flat_alleles())
    print(con.consensus)
    assert len(con.get_flat_alleles()) == 2
    assert con.get_flat_alleles()[0] == "DPA1*01:04"
    assert con.get_flat_alleles()[1] == "DPA1*02"


def test_allele_sorting():
    calls = {
        "sample": "",
        "calls": {
            "C": {
                "alg1": ["C*07:01", "C*12:02"],
                "alg2": ["C*12:02", "C*07:01"],
            }
        },
    }
    report = Report(calls)
    consensus = ConsensusAlgorithm(report)

    alleles = consensus.get_flat_alleles()
    print(alleles)
    assert alleles[0] == "C*07:01"
    assert alleles[1] == "C*12:02"


def test_homozygous_at_lowres():
    calls = make_report([])
    report = Report(calls)

    alleles = [
        "DPA1*02:01:01 (1.0000)",
        "DPA1*02:01:01:02 (0.5329)",
        "DPA1*02:01:01:01 (0.4671)",
    ]

    two_alleles = get_allele_pair(alleles, 2)

    print(two_alleles)

    assert len(two_alleles) == 2
    assert two_alleles[0] == Allele("DPA1*02:01")
    assert two_alleles[1] == Allele("DPA1*02:01")


def test_hetero_at_highres():
    calls = make_report([])
    report = Report(calls)

    alleles = [
        "DPA1*02:01:01 (1.0000)",
        "DPA1*02:01:01:02 (0.5329)",
        "DPA1*02:01:01:01 (0.4671)",
    ]

    two_alleles = get_allele_pair(alleles, 4)

    print(two_alleles)

    assert len(two_alleles) == 2
    assert two_alleles[0] == Allele("DPA1*02:01:01:01")
    assert two_alleles[1] == Allele("DPA1*02:01:01:02")
