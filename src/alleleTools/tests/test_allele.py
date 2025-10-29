import pytest

from typing import List

from ..allele import Allele, AlleleMatchStatus, FieldTree, AlleleParser


def build_allele_tree(gene: str, alleles: List[Allele]) -> FieldTree:
    """
    Build a FieldTree representing the structure of allele fields for a given gene.

    Args:
        gene (str): The gene name to use as the root of the tree.
        alleles (List[Allele]): List of Allele objects to add to the tree.

    Returns:
        FieldTree: The root of the constructed field tree.

    Example:
        >>> build_allele_tree('A', [Allele('A*01:01'), Allele('A*01:02')])
        Field(A:0)[Field(01:2)[Field(01:1), Field(02:1)]]
    """
    root = FieldTree(gene)
    for allele in alleles:
        root.add(allele.fields)
    return root

class TestAlleleParser:
    def test_default_parser(self):
        parser = AlleleParser(gene_family="hla")
        assert "hla" in parser.config
        assert "kir" in parser.config
        assert "hla_hisat" in parser.config
        assert "hla_delimited" in parser.config

    def test_hla_parsing(self):
        parser = AlleleParser(gene_family="hla")
        allele_str = "A*02:01:01:01 (high confidence)"
        allele = parser.parse(allele_str)

        print(allele)
        
        assert allele.gene == "A"
        assert allele.fields == ["02", "01", "01", "01"]

    def test_hla_hisat_parsing(self):
        parser = AlleleParser(gene_family="hla_hisat")
        allele_str = "A*02:01:01:01 (0.5)"
        allele = parser.parse(allele_str)

        print(allele)

        assert allele.gene == "A"
        assert allele.fields == ["02", "01", "01", "01"]
        assert allele.confidence == 0.5

    def test_delimited_parsing(self):
        parser = AlleleParser(gene_family="hla_delimited")
        allele_str = "A*02:01:01:01"
        allele = parser.parse(allele_str)

        print(allele)
        
        assert allele.gene == "A"
        assert allele.fields == ["02", "01", "01", "01"]
    
    def test_wrong_parser(self):
        # check that the exception is raised
        with pytest.raises(Exception):
            parser = AlleleParser(gene_family="wrong_parser")
            allele_str = "InvalidAlleleString"
            parser.parse(allele_str)


class TestHlaAllele:
    parser = AlleleParser(gene_family="hla")

    def test_comparison_results(self):
        a1 = self.parser.parse("A*02:01")
        a2 = self.parser.parse("A*02:01:01:01")

        assert a1.compare(a2) == AlleleMatchStatus.MORE_RESOLUTION
        assert a2.compare(a1) == AlleleMatchStatus.LESS_RESOLUTION
        assert a1.compare(a1) == AlleleMatchStatus.EQUAL

        a3 = self.parser.parse("B*02:01:01:01")
        assert a1.compare(a3) == AlleleMatchStatus.NOT_EQUAL


    def test_allele_parsing(self):
        a1 = self.parser.parse("A*02:01")
        assert str(a1) == "A*02:01"

        a1 = self.parser.parse("A*02")
        assert str(a1) == "A*02"


class TestKirAllele:
    @pytest.mark.parametrize("allele_str", [
        "KIR2DL1*00302",
        "KIR3DS1*013",
        "KIR2DL4*000"
    ])
    def test_allele_parsing(self, allele_str):
        parser = AlleleParser(gene_family="kir")
        a1 = parser.parse(allele_str)

        assert isinstance(a1, Allele)
        assert str(a1) == allele_str


class TestFieldTree:
    def parse_alleles(self, alleles: list) -> list:
        parser = AlleleParser(gene_family="hla")
        return [parser.parse(a) for a in alleles]

    def test_one_allele(self):
        alleles = ["A*01"]
        alleles = self.parse_alleles(alleles)

        tree = build_allele_tree("A", alleles)

        assert tree.support == 1

    def test_1allele_1field(self):
        alleles = ["A*01"]
        alleles = self.parse_alleles(alleles)

        tree = build_allele_tree("A", alleles)
        assert tree.children[0].field == "01"

    def test_1allele_2fields(self):
        alleles = ["A*01:02"]
        alleles = self.parse_alleles(alleles)

        tree = build_allele_tree("A", alleles)
        assert tree.children[0].children[0].field == "02"

    def test_1allele_3fields(self):
        alleles = ["A*01:02:03"]
        alleles = self.parse_alleles(alleles)

        tree = build_allele_tree("A", alleles)
        assert tree.children[0].children[0].children[0].field == "03"

    def test_1allele_4fields(self):
        alleles = ["A*01:02:03:04"]
        alleles = self.parse_alleles(alleles)

        tree = build_allele_tree("A", alleles)

        assert tree.children[0].children[0].children[0].children[0].field == "04"

    def test_allele_counting_fist_level(self):
        alleles = ["A*01:02", "A*01:02"]
        alleles = self.parse_alleles(alleles)
        tree = build_allele_tree("A", alleles)

        assert tree.support == 2

    def test_allele_counting_second_level(self):
        alleles = ["A*01:02", "A*01:02", "A*02:03"]
        alleles = self.parse_alleles(alleles)
        tree = build_allele_tree("A", alleles)

        assert tree.children[0].support == 2

    def test_allele_counting_second_level_minor_allele(self):
        alleles = ["A*01:02", "A*01:02", "A*02:03"]
        alleles = self.parse_alleles(alleles)
        tree = build_allele_tree("A", alleles)

        assert tree.children[1].support == 1

    def test_allele_support_homozygous(self):
        alleles = ['DPA1*01:03:01', 'DPA1*01:03:01:04', 'DPA1*01:03:01:04']
        alleles = self.parse_alleles(alleles)
        tree = build_allele_tree("DPA1", alleles)

        alleles, support = tree.get_consensus(0.6)

        assert alleles == ["DPA1*01:03:01:04", "DPA1*01:03:01:04"]
