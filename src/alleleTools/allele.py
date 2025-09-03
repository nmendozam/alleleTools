import re
from enum import Enum
from math import ceil
from typing import List, Tuple


class AlleleMatchStatus(Enum):
    """
    Enumeration for allele comparison results.

    Defines the possible outcomes when comparing two alleles:
    - NOT_EQUAL: Alleles are different (different genes or field values)
    - EQUAL: Alleles are identical in gene and all fields
    - LESS_RESOLUTION: Current allele has fewer fields than compared allele
    - MORE_RESOLUTION: Current allele has more fields than compared allele
    """

    NOT_EQUAL = 1
    EQUAL = 2
    LESS_RESOLUTION = 3
    MORE_RESOLUTION = 4


class Allele:
    """
    Class representing an HLA allele with parsing and comparison capabilities.

    This class handles the parsing of HLA allele nomenclature and provides
    methods for comparison, resolution management, and string representation.
    It supports confidence scores from algorithms like HiSat.

    Attributes:
        gene (str): Gene name (e.g., 'A', 'B', 'DRB1')
        fields (List[str]): List of allele fields (e.g., ['01', '01', '01'])
        confidence (float): Optional confidence score from genotyping algorithm

    Args:
        code (str): Allele string to parse (e.g., 'A*01:01:01' or 'A*01:01(0.95)')

    Raises:
        Exception: If the allele string cannot be parsed

    Example:
        >>> allele = Allele('A*01:01:01')
        >>> print(allele.gene)  # 'A'
        >>> print(allele.fields)  # ['01', '01', '01']
        >>> print(len(allele))  # 3
    """

    def __init__(self, code: str, gene: str = "") -> None:
        if not code:
            raise Exception("No allele to parse")

        # Extract confidence score (Hisat)
        confidence_score = re.search(r"\((.*?)\)", code)
        if confidence_score:
            self.confidence = float(confidence_score.group(1))

        if gene:
            self.gene = gene
            # Extract allele name 01:01:01
            self.__parse_fields(code)
            return

        # Extract allele name GENE*01:01:01
        self.__parse_gene_and_fields(code)

    def __parse_fields(self, code: str):
        pattern = r"(\d{2})(:\d{2})*"
        allele_pattern = re.search(pattern, code)

        if not allele_pattern:
            raise Exception(f"Error allele '{code}' is not parsable")

        allele_code = allele_pattern.group(0)
        self.fields = allele_code.split(":")

    def __parse_gene_and_fields(self, code: str):
        pattern = r"(\w+)\*(\d{2})(:\d{2})*"
        allele_pattern = re.search(pattern, code)

        if not allele_pattern:
            raise Exception(f"Error allele '{code}' is not parsable")

        allele_code = allele_pattern.group(0)
        self.gene, fields = allele_code.split("*")
        self.fields = fields.split(":")

    def __repr__(self) -> str:
        if not hasattr(self, "gene"):
            return "<Empty Allele>"
        return f"""({len(self.fields)}-fields){str(self)}"""

    def __str__(self) -> str:
        """
        String representation of the allele in standard HLA nomenclature.

        Returns:
            str: Allele in format 'GENE*field1:field2:field3'
        """
        if not hasattr(self, "gene"):
            return ""
        return f"{self.gene}*{':'.join(self.fields)}"

    def __len__(self) -> int:
        """
        Get the resolution level of the allele.

        Returns:
            int: Number of fields in the allele (resolution level)
        """
        return len(self.fields)

    def __eq__(self, allele_b: object) -> bool:
        """Check equality with another allele."""
        if not isinstance(allele_b, Allele):
            return NotImplemented

        return self.compare(allele_b) == AlleleMatchStatus.EQUAL

    def __hash__(self):
        """Hash function for use in sets and dictionaries."""
        return hash(str(self))

    def truncate(self, new_resolution: int):
        """
        Reduce the resolution of the allele to the specified level.

        Args:
            new_resolution (int): Target resolution level (number of fields)

        Note:
            If new_resolution is greater than current resolution, no change is made.
        """
        if new_resolution > len(self):
            return
        self.fields = self.fields[:new_resolution]

    def compare(self, allele: "Allele") -> AlleleMatchStatus:
        """
        Compare this allele with another allele.

        Performs detailed comparison considering gene name, field values,
        and resolution levels.

        Args:
            allele (Allele): The allele to compare with

        Returns:
            CmpResult: The result of the comparison

        Example:
            >>> a1 = Allele('A*01:01')
            >>> a2 = Allele('A*01:01:01')
            >>> a1.compare(a2)  # CmpResult.MORE_RESOLUTION
        """
        # Check that it belongs to the same gene
        if self.gene != allele.gene:
            return AlleleMatchStatus.NOT_EQUAL

        # Compare each field
        for s_field, a_field in zip(self.fields, allele.fields):
            if s_field != a_field:
                return AlleleMatchStatus.NOT_EQUAL

        # Check differences in resolution
        if len(self.fields) > len(allele.fields):
            return AlleleMatchStatus.LESS_RESOLUTION
        elif len(self.fields) < len(allele.fields):
            return AlleleMatchStatus.MORE_RESOLUTION

        return AlleleMatchStatus.EQUAL


class FieldTree:
    """
    Tree structure to represent the hierarchical organization of allele fields.

    Each node in the tree corresponds to a field value at a particular position
    in the allele nomenclature. The tree is used to count and organize the
    occurrence of each field value across a set of alleles.

    Attributes:
        field (str): The field value at this node.
        support (int): Number of times this field value has been added at this position,
            which indicates how many times it has been genotyped by different tools.
        children (List[FieldTree]): List of child nodes representing subsequent fields.

    Example:
        For alleles A*01:01 and A*01:02, the tree will have a root 'A',
        a child '01', and two children '01' and '02' under it.
    """

    def __init__(self, name: str):
        """
        Initialize a FieldTree node.

        Args:
            name (str): The field value for this node.
        """
        self.field = name
        self.support = 0
        self.children = list()

    def __str__(self):
        """
        String representation of the FieldTree node and its children.

        Returns:
            str: Human-readable representation of the node and its subtree.
        """
        if len(self.children) == 0:
            return f"Field({self.field}:{self.support})"

        return f"Field({self.field}:{self.support}){self.children}"

    def __repr__(self):
        """
        Official string representation of the FieldTree node.

        Returns:
            str: String representation (same as __str__).
        """
        return self.__str__()

    def add(self, fields: list):
        """
        Add a sequence of fields to the tree, incrementing counts and creating
        nodes as needed.

        Args:
            fields (list): List of field values (str) to add as a path in the
            tree.

        Example:
            tree.add(['01', '01']) will add/increment nodes for '01' at two
            levels.
        """
        self.support += 1

        if len(fields) == 0:
            return

        name = fields[0]
        overhead = fields[1:]

        for child in self.children:
            if child.field == name:
                child.add(overhead)
                return

        # if nothing was found
        new_tree = FieldTree(name)
        new_tree.add(overhead)
        self.children.append(new_tree)

    def get_consensus(self, min_support: float) -> Tuple[List[str], List[float]]:
        """
        Gets a list of up to two posible consensus solutions that meet the criteria
        of the minimum support. This is basically a tree search algorithm.

        Args:
            min_support (float): minimum proportion of support required. Each program
                contributes votes equally, the amount of times the allele has been genotyped
                is stored in each node of the tree. This value is used to filter consensus
                alleles based on their support values.

        Returns:
            Tuple[List[str], List[float]]: Consensus alleles and their support values.
        """
        # min_support is a proportion between 0 and 1
        assert min_support <= 1 and min_support >= 0

        # Look for consensus solutions
        max_support = self.support
        min_support_n = max_support * min_support
        solutions = self.__get_consensus__(ceil(min_support_n / 2))

        # Get the two most supported alleles
        solutions.sort(key=lambda i: i[1], reverse=True)
        solutions = solutions[:2]

        # Reformat the output
        alleles = [a for a, _ in solutions]
        supports = [s / max_support for _, s in solutions]

        return alleles, supports

    def __get_consensus__(self, min_support_number: int) -> List[Tuple[str, float]]:
        """
        Get consensus solutions from the field tree.

        Args:
            min_support_number (int): Minimum support threshold.

        Returns:
            List[Tuple[str, float]]: List of consensus solutions with their support.
        """
        if self.support < min_support_number:
            return list()
        if len(self.children) <= 0:
            return [(self.field, self.support)]

        solutions = list()
        for child in self.children:
            child_sol = child.__get_consensus__(min_support_number)
            solutions.extend(self.__merge_with_current_node__(child_sol))
        return solutions

    def __merge_with_current_node__(self, res: List[Tuple[str, float]]) -> List[Tuple[str, float]]:
        """
        Takes as input the list of possible consensus from child nodes
        and merges it with the current node's field. If no response was
        obtained from the child nodes, it will return the current node's
        field and its support count.

        example:
            if the current node is "A" and the child node is "B:2", the result will be
            ["A:B:2", 2]
        """
        if len(res) == 0:
            return [(self.field, self.support)]
        sol = list()
        for fields, support in res:
            sol.append((f"{self.field}:{fields}", support))
        return sol


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
