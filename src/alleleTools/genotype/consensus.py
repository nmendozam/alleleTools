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

    # Check if output file exists
    if os.path.exists(args.out):
        print(f"Output file {args.out} already exists.")
        exit(1)

    for file in args.input:
        report, consensus = _open_report(file)
        alleles = [str(c) for gene in consensus for c in consensus[gene].alleles]

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
            with open(args.out, "a") as out:
                out.write("%s\t2\t%s\n" % (report.sample, "\t".join(alleles)))

    if args.format == "vcf":
        with open(args.out, "w") as vcf:
            vcf.write(VCF_HEADER)

        df.fillna("0|0:0:1,0,0", inplace=True)
        df.sort_index(inplace=True)
        df.to_csv(args.out, sep="\t", mode="a", index=False)


class ComparisonResult(Enum):
    """Codes for comparison of alleles"""

    NOT_EQUAL = 1
    EQUAL = 2
    LESS_RESOLUTION = 3
    MORE_RESOLUTION = 4


class AlleleParser:
    """Pareses and compares alleles"""

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
        return f"{self.gene}*{':'.join(self.fields)}"


class ReportParser:
    def __init__(self, report) -> None:
        calls = report["calls"]
        self.sample = report["sample"]
        self.genes = dict()
        for gene, call in calls.items():
            self.genes[gene] = self.parse_call(call)

    def parse_call(self, call):
        calls = list()
        flat_list = [
            (algorithm, called_alleles)
            for algorithm in call
            for called_alleles in call[algorithm]
        ]
        for algorithm, pred in flat_list:
            try:
                allele = AlleleParser(algorithm, pred)
                calls.append(allele)
            except:
                pass
        return calls


class ConsensusAlgorithm:
    def __init__(self, gene_name: str) -> None:
        self.gene_name = gene_name
        self.alleles = list()

    def __repr__(self) -> str:
        return f"{self.gene_name}->{self.alleles}"

    class Consensus:
        def __init__(self, new) -> None:
            self.gene = new.gene
            self.fields = new.fields
            self.evidence = [new]

        def __repr__(self) -> str:
            return f"{self.gene}*{':'.join(self.fields)}"
            # return f"{self.gene}*{':'.join(self.fields)} ({len(self.evidence)})"

        def is_equal(self, allele: AlleleParser):
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

    def add(self, new_pred):
        no_matches = True
        for pred in self.alleles:
            result = pred.is_equal(new_pred)
            if result == ComparisonResult.NOT_EQUAL:
                continue

            pred.evidence.append(new_pred)
            no_matches = False

            if result == ComparisonResult.MORE_RESOLUTION:
                pred.fields = new_pred.fields

        if no_matches:
            self.alleles.append(self.Consensus(new_pred))


def _open_report(file_name: str):
    with open(file_name) as f:
        json_report = json.load(f)

        print(file_name)
        report = ReportParser(json_report)
        consensus = dict()
        for gene, calls in report.genes.items():
            consensus[gene] = ConsensusAlgorithm(gene)
            for call in calls:
                consensus[gene].add(call)

            # Print consensus that have more than one allele
            if len(consensus[gene].alleles) > 2:
                for call in calls:
                    print("\t", call, call.program_name)
                print("consensus:", consensus[gene])
                consensus[gene].alleles = consensus[gene].alleles[:2]
            elif len(consensus[gene].alleles) == 1:
                consensus[gene].alleles = consensus[gene].alleles * 2
            elif len(consensus[gene].alleles) == 0:
                consensus[gene].alleles = ["NA", "NA"]

            assert len(consensus[gene].alleles) == 2
        return report, consensus
