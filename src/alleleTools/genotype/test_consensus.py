from alleleTools.allele import AlleleParser
from mock import patch

from .consensus import ConsensusGene


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


def test_consensus():
    with patch.object(ConsensusGene, "__init__", lambda self, *args: None):
        con_gene = ConsensusGene.__new__(ConsensusGene)
        con_gene.allele_parser = AlleleParser(gene_family="hla")
        con_gene.name = "DRB1"
        con_gene.calls = {
            "HLA-HD": ["15:01:01", "04:01:01"],
            "xHLA": ["04:01", "15:01"],
            "HLAscan": ["04:01:01", "15:01:01:04"],
            "Optitype": [],
            "Hisat": [
                "15:01:01 (0.6466)",
                "04:01:01 (0.3534)",
                "15:01:01:02 (0.3443)",
                "15:01:01:01 (0.3023)",
            ],
        }
        alleles, support = con_gene.get_consensus_call(min_support=0.6)
        assert len(alleles) == 2
        assert alleles[0] == "DRB1*04:01:01"
        assert support[0] == 3
        assert alleles[1] == "DRB1*15:01:01"
        assert support[1] == 3


def test_overlapping_allele_consensus():
    with patch.object(ConsensusGene, "__init__", lambda self, *args: None):
        con_gene = ConsensusGene.__new__(ConsensusGene)
        con_gene.allele_parser = AlleleParser(gene_family="hla")
        con_gene.name = "DPA1"
        con_gene.calls = {
            "alg1": [
                "DPA1*01 (0.6666)",
                "DPA1*01:03:01 (0.3333)",
                "DPA1*01:04 (0.3333)",
                "DPA1*01:03:01:01 (0.3114)",
            ]
        }
        alleles, support = con_gene.get_consensus_call(min_support=0.6)
        assert len(alleles) == 2
        assert alleles[0] == "DPA1*01:03:01:01"
        assert alleles[1] == "DPA1*01:04"


def test_allele_sorting():
    with patch.object(ConsensusGene, "__init__", lambda self, *args: None):
        con_gene = ConsensusGene.__new__(ConsensusGene)
        con_gene.name = "C"
        con_gene.calls = {
            "alg1": ["C*07:01", "C*12:02"],
            "alg2": ["C*12:02", "C*07:01"],
        }
        con_gene.allele_parser = AlleleParser(gene_family="hla")
        alleles, support = con_gene.get_consensus_call(min_support=0.6)
        assert alleles[0] == "C*07:01"
        assert alleles[1] == "C*12:02"