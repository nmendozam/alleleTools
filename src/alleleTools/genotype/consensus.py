import argparse
import glob
import json
import os
import re
from enum import Enum

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

        consensus = ConsensusAlgorithm(report)
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


class ComparisonResult(Enum):
    """Codes for comparison of alleles"""

    NOT_EQUAL = 1
    EQUAL = 2
    LESS_RESOLUTION = 3
    MORE_RESOLUTION = 4


class Allele:
    """Class to represent alleles"""

    def __init__(self, program_name: str, code: str) -> None:
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
        self.program_name = program_name

    def __repr__(self) -> str:
        if not hasattr(self, "gene"):
            return ""
        # return f"({len(self.fields)}-fields){self.gene}*{':'.join(self.fields)}"
        return f"{self.gene}*{':'.join(self.fields)}"

    def compare(self, allele: "Allele") -> ComparisonResult:
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
    def __init__(self, report: dict) -> None:
        calls = report["calls"]
        self.sample = report["sample"]
        self.genes = dict()
        for gene, call in calls.items():
            self.genes[gene] = self.parse_call(call)

    def parse_call(self, call: dict) -> list[Allele]:
        alleles = list()

        for algorithm, called_alleles in call.items():
            for allele in called_alleles:
                try:
                    allele = Allele(algorithm, allele)
                    alleles.append(allele)
                except Exception as e:
                    print(f"Error parsing allele {allele} from {algorithm} algorithm")
                    print("\t%s" % e)
                    pass
        return alleles


class ConsensusAlgorithm:
    def __init__(self, report: Report) -> None:
        self.consensus = dict()

        for gene, calls in report.genes.items():
            alleles = list()

            # Cluster similar alleles
            for call in calls:
                match = self.find_matching_allele(call, alleles)

                if not match:
                    alleles.append(self.Consensus(call))
                else:
                    match.add_evidence(call)

                    if match.compare(call) == ComparisonResult.MORE_RESOLUTION:
                        match.fields = call.fields

            self.consensus[gene] = alleles

    def correct_homozygous_calls(self):
        """Go over alleles and correct homozygosity"""
        for gene in self.consensus:
            n_alleles = len(self.consensus[gene])

            if n_alleles > 2:
                self.consensus[gene] = self.consensus[gene][:2]
            elif n_alleles == 1:
                self.consensus[gene] = self.consensus[gene] * 2
            elif n_alleles == 0:
                self.consensus[gene] = ["NA", "NA"]

    def find_matching_allele(self, allele: Allele, allele_list):
        """Finds a matching allele in the consensus"""
        for a in allele_list:
            result = a.compare(allele)
            if result != ComparisonResult.NOT_EQUAL:
                return a
        return None

    def get_flat_alleles(self) -> list[str]:
        return [str(c) for gene in self.consensus for c in self.consensus[gene]]

    class Consensus:
        def __init__(self, new) -> None:
            self.main_allele: Allele = new
            self.evidence = [new]

        def __repr__(self) -> str:
            return str(self.main_allele)
            # return f"{self.main_allele} ({len(self.evidence)})"
        
        def add_evidence(self, allele: Allele) -> None:
            self.evidence.append(allele)
        
        def compare(self, allele: Allele) -> ComparisonResult:
            return self.main_allele.compare(allele)


