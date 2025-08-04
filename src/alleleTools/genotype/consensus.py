import argparse
import json
import re
from enum import Enum
from typing import Dict, List

import pandas as pd

VCF_ROW_REF = {
    "CHROM": [6],
    "POS": [],
    "ID": [],
    "REF": ["A"],
    "ALT": ["T"],
    "QUAL": ["."],
    "FILTER": ["PASS"],
    "INFO": ["AR2=1.00;DR2=1.00;AF=0.16"],
    "FORMAT": ["GT:DS:GP"],
}

VCF_HEADER = """##fileformat=VCFv4.1
##fileDate=20090805
##source=alleleTools_consensus
##reference=file:///seq/references/
#"""


def setup_parser(subparsers):
    parser = subparsers.add_parser(
        name="consensus",
        description="This program finds a consensus between multiple HLA genotyping reports",
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
        type=str,
        help="Path to output file",
        default="output.alt",
    )
    parser.add_argument(
        "--format",
        choices=["vcf", "alt"],
        default="alt",
        type=str,
        help="Format of the output file",
    )

    parser.set_defaults(func=call_function)

    return parser


def call_function(args):
    if args.format == "vcf":
        # Relevant variables for vcf format
        df = pd.DataFrame(columns=VCF_ROW_REF.keys())
        position = 29910247

    for file in args.input:
        json_report = json.load(file)

        report = Report(json_report)

        consensus = ConsensusAlgorithm(report.genes)
        consensus.correct_homozygous_calls()
        alleles = consensus.get_flat_alleles()

        ##########################
        # Generate the output file
        ##########################

        if args.format == "vcf":
            # Look for existing alleles
            for allele in alleles[:]:
                if allele in df.index:
                    df.loc[allele, report.sample] = "1|0:1:0,1,0"
                    alleles.remove(allele)

            # Add new alleles
            vcf_row = {k: v * len(alleles) for k, v in VCF_ROW_REF.items()}
            vcf_row["POS"] = range(position, position + len(alleles))
            vcf_row["ID"] = alleles
            vcf_row[report.sample] = ["1|0:1:0,1,0"] * len(alleles)

            position += len(alleles)
            # Add consensus to data frame
            df = pd.concat([df, pd.DataFrame(vcf_row, index=alleles)])
        elif args.format == "alt":
            with open(args.output, "a") as out:
                out.write("%s\t2\t%s\n" % (report.sample, "\t".join(alleles)))

    if args.format == "vcf":
        with open(args.output, "w") as vcf:
            vcf.write(VCF_HEADER)

        df.fillna("0|0:0:1,0,0", inplace=True)
        df.sort_index(inplace=True)
        df.to_csv(args.output, sep="\t", mode="a", index=False)


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
            abundance += allele.confidence if hasattr(allele, "confidence") else 0.5
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


class ComparisonResult(Enum):
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
            return ""
        return f"({len(self.fields)}-fields){self.gene}*{':'.join(self.fields)}"

    def __str__(self) -> str:
        """String representation of the allele"""
        if not hasattr(self, "gene"):
            return ""
        return f"{self.gene}*{':'.join(self.fields)}"

    def __len__(self) -> int:
        """Returns the resolution of the allele in n-fields"""
        return len(self.fields)

    def __eq__(self, allele_b: "Allele"):
        return self.compare(allele_b) == ComparisonResult.EQUAL
    
    def __hash__(self):
        return hash(str(self))

    def truncate(self, new_resolution: int):
        """Reduce the resolution of the allele"""
        if new_resolution > len(self):
            return
        self.fields = self.fields[:new_resolution]

    def compare(self, allele: "Allele") -> ComparisonResult:
        """Compares this allele to another. Check if they are equal, not equal, or if one has more resolution than the other."""
        # Check that it belongs to the same gene
        if self.gene != allele.gene:
            return ComparisonResult.NOT_EQUAL

        # Compare each field
        for s_field, a_field in zip(self.fields, allele.fields):
            if s_field != a_field:
                return ComparisonResult.NOT_EQUAL

        # Check differences in resolution
        if len(self.fields) > len(allele.fields):
            return ComparisonResult.LESS_RESOLUTION
        elif len(self.fields) < len(allele.fields):
            return ComparisonResult.MORE_RESOLUTION

        return ComparisonResult.EQUAL


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


class ConsensusAlgorithm:
    def __init__(self, calls: Dict[str, Dict[str, List[Allele]]]) -> None:
        self.consensus = dict()

        for gene, alleles in calls.items():
            allele_cluster = list()

            # Cluster similar alleles
            for call in alleles.values():
                for allele in call:
                    match = self.find_matching_allele(allele, allele_cluster)

                    if not match:
                        allele_cluster.append(self.AlleleWithEvidence(allele))
                    else:
                        match.add_evidence(allele)

                        if match.compare(allele) == ComparisonResult.MORE_RESOLUTION:
                            match.fields = allele.fields

            self.consensus[gene] = allele_cluster

    def correct_homozygous_calls(self):
        """Go over alleles and correct zygosity. e.i. if there is only one allele,
        the allele will be duplicated. If there is more than two, the allele list
        will be truncated. If the list is empty, it will return NA"""
        for gene in self.consensus:
            n_alleles = len(self.consensus[gene])

            if n_alleles > 2:
                self.consensus[gene] = self.consensus[gene][:2]
            elif n_alleles == 1:
                self.consensus[gene] = self.consensus[gene] * 2
            elif n_alleles == 0:
                self.consensus[gene] = ["NA", "NA"]

    def find_matching_allele(
        self, allele: Allele, allele_list: list["AlleleWithEvidence"]
    ) -> "AlleleWithEvidence":
        """Finds a matching allele in the consensus"""
        for a in allele_list:
            result = a.compare(allele)
            if result != ComparisonResult.NOT_EQUAL:
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

