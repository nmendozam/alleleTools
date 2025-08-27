import re
from enum import Enum


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

    def __init__(self, code: str) -> None:
        if not code:
            raise Exception("No allele to parse")

        # Extract allele name GENE*01:01:01
        allele_pattern = re.search(r"(\w+)\*(\d{2})(:\d{2})*", code)
        if not allele_pattern:
            raise Exception("Error no allele name parsable")

        # Extract confidence score (Hisat)
        confidence_score = re.search(r"\((.*?)\)", code)
        if confidence_score:
            self.confidence = float(confidence_score.group(1))

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

    def __eq__(self, allele_b: "Allele"):
        """Check equality with another allele."""
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

    def compare(self, allele: "Allele") -> CmpResult:
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
        return AlleleMatchStatus.EQUAL

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