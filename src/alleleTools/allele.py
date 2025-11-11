import json
import math
import re
from enum import Enum
from typing import Iterator, List, Tuple, Union
from abc import ABC, abstractmethod

from .utils.assets import get_asset_path


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
        code (str | list): Allele string to parse (e.g., 'A*01:01:01' or
            'A*01:01(0.95)'). Alternatively, you can pass a list of fields
            along with the argument `gene`
        gene (str): If the gene name is not present in `code` you should pass
            the here here.

    Raises:
        Exception: If the allele string cannot be parsed

    Example:
        >>> allele = Allele('A*01:01:01')
        >>> print(allele.gene)  # 'A'
        >>> print(allele.fields)  # ['01', '01', '01']
        >>> print(len(allele))  # 3
    """

    def __init__(
        self,
        gene: str,
        fields: List[str],
        confidence: str | None = None,
        gene_delimiter: str = "*",
        field_delimiter: str = ":",
    ):
        self.fields: List[str] = fields
        self.gene = gene
        self.gene_delimiter = gene_delimiter
        self.field_delimiter = field_delimiter

        # Extract confidence score (Hisat)
        if confidence:
            self.confidence = float(confidence)

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
        if not self.gene and not self.fields:
            return ""
        return (
            f"{self.gene}{self.gene_delimiter}{self.field_delimiter.join(self.fields)}"
        )

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

    def get_fields(self) -> List[str]:
        return self.fields

    def truncate(self, new_resolution: int) -> "Allele":
        """
        Reduce the resolution of the allele to the specified level.

        Args:
            new_resolution (int): Target resolution level (number of fields)

        Note:
            If new_resolution is greater than current resolution, no change is made.
        """
        if new_resolution > len(self):
            return Allele("", [])
        self.fields = self.fields[:new_resolution]
        return self

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


class ParsingStrategy(ABC):
    """
    Abstract base class for allele parsing strategies.
    """

    @abstractmethod
    def parse(self, text: str) -> Allele:
        pass


class DelimitedParser(ParsingStrategy):
    """
    Parser for alleles using simple delimiters. Two parameters are needed:
        - gene_delimiter: character that separates the gene name from the fields
        - field_delimiter: character that separates the different fields

    Example:
        For allele "A*01:02:03", gene_delimiter="*", field_delimiter=":"
    Results:
        1. gene = "A"
        2. fields = ["01", "02", "03"]
    """

    def __init__(self, gene_delimiter: str, field_delimiter: str):
        self.gene_delimiter = gene_delimiter
        self.field_delimiter = field_delimiter

    def parse(self, text: str) -> Allele:
        parts = text.split(self.gene_delimiter)
        if len(parts) != 2:
            print(f"Warning: Could not parse allele from text '{text}'")
            return Allele("", [])
        gene = parts[0]
        field_part = parts[1]
        fields = field_part.split(self.field_delimiter)
        return Allele(
            gene=gene,
            fields=fields,
            gene_delimiter=self.gene_delimiter,
            field_delimiter=self.field_delimiter,
        )


class RegexParser(ParsingStrategy):
    """
    Parser for alleles using regular expressions. It only requires the regex pattern.
    It should contain named groups for 'gene', 'field1', 'field2', etc. An optional
    group 'confidence' can also be included to capture confidence scores.

    field_delimiter and gene_delimiter can also be specified. These parameters will only be used
    to format the the allele string representation.

    This parser is more flexible and robust than the DelimitedParser, as it can handle
    more complex allele formats. However, it requires knowledge of regular expressions.

    Example:
        For allele "A*01",
        pattern = r"(?:(?P<gene>\\w+)\\*)?(?P<field1>\\d{2})"
    Results:
        1. gene = "A"
        2. fields = ["01"]
    """

    def __init__(self, pattern: str, field_delimiter: str, gene_delimiter: str):
        self.pattern = re.compile(pattern)
        self.field_delimiter = field_delimiter
        self.gene_delimiter = gene_delimiter

    def parse(self, text: str) -> Allele:
        match = self.pattern.search(text)

        if not match:
            print(f"Warning: Could not parse allele from text '{text}'")
            return Allele("", [])

        fields = []
        result = match.groupdict()
        for key, value in result.items():
            if key.startswith("field") and value is not None:
                fields.append(value)

        return Allele(
            gene=result.get("gene", ""),
            fields=fields,
            confidence=result.get("confidence", None),
            field_delimiter=self.field_delimiter,
            gene_delimiter=self.gene_delimiter,
        )


class AlleleParser:
    """
    Configurable allele parser that supports multiple parsing strategies.
    This class loads parsing configurations from a JSON file.

    The configuration can be overwritten by providing a custom config file during
    initialization (config_file parameter).

    Use:
        parser = AlleleParser(gene_family="hla", config_file="custom_config.json")
        allele = parser.parse("A*01:02:03")
    """

    def __init__(self, gene_family: str, config_file: str = ""):
        # Load default config
        config_def = get_asset_path("parser_config.json")
        with open(config_def, "r") as f:
            self.config = json.load(f)

        # Load custom config
        if config_file != "":
            with open(config_file, "r") as f:
                self.config.update(json.load(f))

        self.strategies = self._build_strategies()
        self.set_gene_family(gene_family)

    def _build_strategies(self):
        strategies = dict()
        for name, config in self.config.items():
            if config["type"] == "delimited":
                strategies[name] = DelimitedParser(
                    config["gene_delimiter"], config["field_delimiter"]
                )
            elif config["type"] == "regex":
                strategies[name] = RegexParser(
                    config["pattern"],
                    field_delimiter=config.get("field_delimiter", ":"),
                    gene_delimiter=config.get("gene_delimiter", "*"),
                )

        return strategies

    def set_gene_family(self, gene_family: str):
        if gene_family not in self.strategies:
            raise Exception(
                f"Gene '{gene_family}' not found in allele parser configuration"
            )
        self.gene_family = gene_family

    def get_delimiters(self) -> Tuple[str, str]:
        strategy = self.strategies[self.gene_family]
        return strategy.gene_delimiter, strategy.field_delimiter

    def parse(self, code: str) -> Allele:
        return self.strategies[self.gene_family].parse(code)


class Solution:
    def __init__(self, allele: str, support: float) -> None:
        self.allele = allele
        self.support = support


class FieldTree:
    """
    Tree structure to represent the hierarchical organization of allele fields.

    Each node in the tree corresponds to a field value at a particular position
    in the allele nomenclature. The tree is used to count and organize the
    occurrence of each field value across a set of alleles.

    Attributes:
        field (str): The field value at this node.
        support (int): Number of times this field value has been added at this
            position, which indicates how many times it has been genotyped by
            different tools.
        children (List[FieldTree]): List of child nodes
        representing subsequent fields.

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

    def set_support(self, new_weight: float, recursive: bool = True) -> None:
        """
        Change the support value. When recursive is True, it will change the
        support value for all children.
        """
        self.support = new_weight

        if not recursive:
            return

        for child in self.children:
            child.set_support(new_weight, recursive=recursive)

    def add_batch(self, batch: List[list], weight: float = 1.0):
        """
        Helper method to append a list of field lists to the current tree.

        Args:
            batch (List[list]): A list containing lists of field values (str).
        """
        if len(batch) == 0:
            return

        for fields in batch:
            self.add(fields, weight=weight)

    def merge_tree(self, tree: "FieldTree"):
        """
        Merges a foreign tree into the current tree. The top level should
        match with this tree.

        Args:
            tree (FieldTree): The foreign tree.
        """
        assert tree.field == self.field

        self.support += tree.support

        children = {child.field: child for child in self.children}

        for t_child in tree.children:
            if t_child.field not in children.keys():
                self.children.append(t_child)
                continue

            children[t_child.field].merge_tree(t_child)

    def add(self, fields: list, weight: float = 1.0):
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
        self.support += weight

        if len(fields) == 0:
            return

        name = fields[0]
        overhead = fields[1:]

        for child in self.children:
            if child.field == name:
                child.add(overhead, weight)
                return

        # if nothing was found
        new_tree = FieldTree(name)
        new_tree.add(overhead, weight)
        self.children.append(new_tree)

    def get_consensus(
        self,
        min_support: float,
        gene_delimiter: str = "*",
        field_delimiter: str = ":",
        max_support: float = 0,
    ) -> Tuple[List[str], List[float]]:
        """
        Gets a list of up to two possible consensus solutions that meet the
        criteria of the minimum support. This is basically a tree search
        algorithm.

        Two additional parameters can be provided to format the allele strings.

        Args:
            min_support (float): minimum proportion of support required. Each
                program contributes votes equally, the amount of times the
                allele has been genotyped is stored in each node of the tree.
                This value is used to filter consensus alleles based on their
                support values.
            gene_delimiter (str): Delimiter between gene name and fields.
            field_delimiter (str): Delimiter between fields.

        Returns:
            Tuple[List[str], List[float]]: Consensus alleles and their support
                values.
        """
        # min_support is a proportion between 0 and 1
        assert min_support <= 1 and min_support >= 0

        # Look for consensus solutions
        if not max_support:
            max_support = self.support
        min_support_n = max_support * min_support
        solutions = self.__get_consensus__(math.ceil(min_support_n))

        # Get the two most supported alleles
        solutions.sort(key=lambda sol: sol.support, reverse=True)
        solutions = solutions[:2]

        # Now sort by allele name
        solutions.sort(key=lambda sol: sol.allele)

        # For homozygous calls adjust the minimum support threshold
        if len(solutions) == 1 and solutions[0].support >= max_support / 2:
            solutions *= 2

        return self.__format_solutions_as_alleles__(
            solutions, gene_delimiter, field_delimiter
        )

    def __format_solutions_as_alleles__(
        self, solutions: List[Solution], gene_delimiter: str, field_delimiter: str
    ) -> Tuple[List[str], List[float]]:
        """
        Formats the solutions as two separate lists:

        Args:
            solutions (list): is a list of the found consensus alleles with
                their scores
            gene_delimiter (str): Delimiter between gene name and fields.
            field_delimiter (str): Delimiter between fields.


        Returns:
            - list of consensus alleles (str)
            - list of respective number of algorithms supporting the consensus
              (float)
        """
        alleles, supports = list(), list()
        for solution in solutions:
            supports.append(solution.support)

            nodes = solution.allele.split(":")

            # If it's only the gene name, there was no consensus
            if len(nodes) <= 1:
                alleles.append("")
                continue

            # Use the allele class as an interface to format the string
            alleles.append(
                str(
                    Allele(
                        gene=nodes[0],
                        fields=nodes[1:],
                        gene_delimiter=gene_delimiter,
                        field_delimiter=field_delimiter,
                    )
                )
            )

        return alleles, supports

    def __get_consensus__(self, min_support_number: float) -> List[Solution]:
        """
        Get consensus solutions from the field tree.

        Args:
            min_support_number (float): Minimum support threshold.

        Returns:
            List[Tuple[str, float]]: List of consensus solutions with their
                support.
        """
        if self.support < min_support_number:
            return list()
        if len(self.children) <= 0:
            return [Solution(self.field, self.support)]

        solutions = list()
        for child in self.children:
            child_sol = child.__get_consensus__(min_support_number)
            if not child_sol:
                continue
            solutions.extend(self.__merge_with_current_node__(child_sol))

        if solutions:
            return solutions

        return [Solution(self.field, self.support)]

    def __merge_with_current_node__(
        self, res: List[Solution]
    ) -> List[Solution]:
        """
        Takes as input the list of possible consensus from child nodes
        and merges it with the current node's field. If no response was
        obtained from the child nodes, it will return the current node's
        field and its support count.

        example:
            if the current node is "A" and the child node is "B:2", the result
            will be ["A:B:2", 2]
        """
        if len(res) == 0:
            return [Solution(self.field, self.support)]
        sol = list()
        for solution in res:
            fields = solution.allele
            support = solution.support
            sol.append(Solution(f"{self.field}:{fields}", support))
        return sol
