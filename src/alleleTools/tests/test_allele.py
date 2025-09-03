from typing import List

from ..allele import Allele, AlleleMatchStatus, FieldTree


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

class TestAllele:
    def test_comparison_results(self):
        a1 = Allele("A*02:01")
        a2 = Allele("A*02:01:01:01")

        assert a1.compare(a2) == AlleleMatchStatus.MORE_RESOLUTION
        assert a2.compare(a1) == AlleleMatchStatus.LESS_RESOLUTION
        assert a1.compare(a1) == AlleleMatchStatus.EQUAL

        a3 = Allele("B*02:01:01:01")
        assert a1.compare(a3) == AlleleMatchStatus.NOT_EQUAL


    def test_allele_parsing(self):
        a1 = Allele("A*02:01")
        assert str(a1) == "A*02:01"

        a1 = Allele("A*02")
        assert str(a1) == "A*02"


    def test_allele_parsing_no_allele(self):
        failed = False
        try:
            a1 = Allele("")
        except:
            failed = True
        assert failed


class TestFieldTree:
    def parse_alleles(self, alleles: list) -> list:
        return [Allele(a) for a in alleles]

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

        assert tree.children[0].num == 2

    def test_allele_counting_second_level_minor_allele(self):
        alleles = ["A*01:02", "A*01:02", "A*02:03"]
        alleles = self.parse_alleles(alleles)
        tree = build_allele_tree("A", alleles)

        assert tree.children[1].num == 1

    def test_allele_support_homozygous(self):
        alleles = ['DPA1*01:03:01', 'DPA1*01:03:01:04', 'DPA1*01:03:01:04']
        alleles = self.parse_alleles(alleles)
        tree = build_allele_tree("DPA1", alleles)

        alleles, support = tree.get_consensus(0.6)

        assert tree.children[1].num == 1

        assert alleles == ["DPA1*01:03:01", "DPA1*01:03:01:04"]
        assert support == [1.0, 0.6666666666666666]
