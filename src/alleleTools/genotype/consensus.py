import argparse
import json
import re
from collections import defaultdict
from enum import Enum
from typing import Dict, List

from ..argtypes import output_path


def setup_parser(subparsers):
    parser = subparsers.add_parser(
        name="consensus",
        description="""
        This program finds a consensus between multiple HLA genotyping reports
        """,
        epilog="Author: Nicolás Mendoza Mejía (2023)",
    )
    parser.add_argument(
        "input",
        metavar="path",
        type=argparse.FileType("r"),
        nargs="+",
        help="JSON files with HLA genotyping reports to be processed",
    )
    parser.add_argument(
        "--output",
        metavar="path",
        type=output_path,
        help="Path to output file",
        default="output.alt",
    )

    parser.set_defaults(func=call_function)

    return parser


def call_function(args):
    for file in args.input:
        json_report = json.load(file)

        report = Report(json_report, resolution=2)

        consensus = ConsensusAlgorithm(report)
        alleles = consensus.get_flat_alleles()

        ##########################
        # Generate the output file
        ##########################

        with open(args.output, "a") as out:
            out.write("%s\t2\t%s\n" % (report.sample, "\t".join(alleles)))


def pairwise(iterable):
    """
    Iterate over a list in a pairwise manner:
        pairwise('ABCDEFG') → AB BC CD DE EF FG
    """

    iterator = iter(iterable)
    a = next(iterator, None)

    for b in iterator:
        yield a, b
        a = b


def get_allele_pair(alleles: List[str], resolution) -> List["Allele"]:
    """Gets a list of alleles in str and returns a list of two
    Allele objects with the desired resolution"""
    parsed_alleles = list()
    alleles.sort()
    for allele in alleles:
        try:
            parsed_allele = Allele(allele)
            parsed_alleles.append(parsed_allele)
        except Exception as e:
            print(f"Error: exeption {e}")
            print(f"The allele {allele} was not parsed")
            pass

    ##############################################################
    # Hisat specific section:
    # Given that Hisat returns a list of most abundant alleles,
    # this sections trys to use the ranking to come up with the
    # highest resolution and most abundant alleles
    ##############################################################

    # Check the abundance of alleles and discard if more than two
    pair = list()
    for pair in pairwise(parsed_alleles):
        abundance = 0
        for allele in pair:
            abundance += allele.confidence if hasattr(
                allele, "confidence") else 0.5
        resolution_is_met = all([len(a) > resolution for a in pair])

        if abundance > 0.90 and resolution_is_met:
            break

    ##############################################################
    # End of Hisat specific section
    ##############################################################

    # truncate the alleles to the desired resolution
    for allele in pair:
        allele.truncate(resolution)

    assert len(pair) <= 2, "More than two alleles found in call"

    return pair[:2]  # Return only the first two alleles, if any


class CmpResult(Enum):
    """Codes for comparison of alleles"""

    NOT_EQUAL = 1
    EQUAL = 2
    LESS_RESOLUTION = 3
    MORE_RESOLUTION = 4


class Allele:
    """Class to represent alleles"""

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
        """String representation of the allele"""
        if not hasattr(self, "gene"):
            return ""
        return f"{self.gene}*{':'.join(self.fields)}"

    def __len__(self) -> int:
        """Returns the resolution of the allele in n-fields"""
        return len(self.fields)

    def __eq__(self, allele_b: "Allele"):
        return self.compare(allele_b) == CmpResult.EQUAL

    def __hash__(self):
        return hash(str(self))

    def truncate(self, new_resolution: int):
        """Reduce the resolution of the allele"""
        if new_resolution > len(self):
            return
        self.fields = self.fields[:new_resolution]

    def compare(self, allele: "Allele") -> CmpResult:
        """
        Compares this allele to another. Check if they are equal, not equal,
        or if one has more resolution than the other.
        """
        # Check that it belongs to the same gene
        if self.gene != allele.gene:
            return CmpResult.NOT_EQUAL

        # Compare each field
        for s_field, a_field in zip(self.fields, allele.fields):
            if s_field != a_field:
                return CmpResult.NOT_EQUAL

        # Check differences in resolution
        if len(self.fields) > len(allele.fields):
            return CmpResult.LESS_RESOLUTION
        elif len(self.fields) < len(allele.fields):
            return CmpResult.MORE_RESOLUTION

        return CmpResult.EQUAL


class Report:
    def __init__(self, report: dict, resolution: int = 2) -> None:
        self.resolution = resolution  # desired resolution in fields
        calls = report["calls"]
        self.sample = report["sample"]
        self.genes = dict()

        for gene, call in calls.items():
            self.genes[gene] = self.parse_call(call)

    def parse_call(self, call: dict) -> dict[Allele]:
        alleles_by_alg = dict()

        for algorithm, alleles in call.items():
            two_alleles = get_allele_pair(alleles, self.resolution)
            alleles_by_alg[algorithm] = two_alleles
        return alleles_by_alg

    def __iter__(self):
        for gene, prog_calls in self.genes.items():
            for program, calls in prog_calls.items():
                for allele in calls:
                    yield {
                        "gene": gene,
                        "program": program,
                        "allele": allele
                    }


class ConsensusAlgorithm:
    def __init__(self, calls: Report) -> None:
        self.consensus = defaultdict(list)

        for report in calls:
            allele = report["allele"]
            gene = report["gene"]
            allele_cluster = self.consensus[gene]
            match = self.find_matching_allele(allele, allele_cluster)

            if not match:
                allele_cluster.append(self.AlleleWithEvidence(allele))
            else:
                match.add_evidence(allele)

                if match.compare(allele) == CmpResult.MORE_RESOLUTION:
                    match.fields = allele.fields

        self.correct_homozygous_calls(self.consensus)

    def correct_homozygous_calls(self, consensus: dict):
        """
        Go over alleles and correct zygosity. e.i. if there is only one
        allele, the allele will be duplicated. If there is more than two,
        the allele list will be truncated. If the list is empty, it will
        return NA
        """
        for gene in consensus:
            n_alleles = len(consensus[gene])

            if n_alleles > 2:
                consensus[gene] = consensus[gene][:2]
            elif n_alleles == 1:
                consensus[gene] = consensus[gene] * 2
            elif n_alleles == 0:
                consensus[gene] = ["NA", "NA"]

    def find_matching_allele(
        self, allele: Allele, allele_list: list["AlleleWithEvidence"]
    ) -> "AlleleWithEvidence":
        """Finds a matching allele in the consensus"""
        for a in allele_list:
            result = a.compare(allele)
            if result != CmpResult.NOT_EQUAL:
                return a
        return None

    def get_flat_alleles(self) -> list[str]:
        return [str(c) for gene in self.consensus for c in self.consensus[gene]]

    class AlleleWithEvidence(Allele):
        def __init__(self, allele: Allele) -> None:
            # Make a copy of the allele's attributes
            self.gene = allele.gene
            self.fields = allele.fields
            if hasattr(allele, "confidence"):
                self.confidence = allele.confidence

            self.evidence = [allele]

        def add_evidence(self, allele: Allele) -> None:
            self.evidence.append(allele)
