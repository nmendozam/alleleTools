"""
kir-mapper to allele table conversion module.

This module reads the reports of kir-mapper and generates
the allele table. Some filtering based on depth and allele
mismatch can be performed.

Author: Nicolás Mendoza Mejía (2025)
"""

from typing import List

import pandas as pd

from alleleTools.allele import AlleleParser
from alleleTools.format.alleleTable import AlleleTable
from alleleTools.format.from_ikmb_hla import ConsensusGene

from ..argtypes import add_out_altable_args, csv_file, file_path, output_path


def setup_parser(subparsers):
    """
    Set up the argument parser for the kir-mapper command.

    Args:
        subparsers: The subparsers object to add this command to.

    Returns:
        argparse.ArgumentParser: The configured parser for kir-mapper.
    """
    parser = subparsers.add_parser(
        name="from_kirmapper",
        help="Convert kir-mapper reports to allele table format",
        description="Convert kir-mapper reports to allele table format",
        epilog="Author: Nicolás Mendoza Mejía (2025)",
    )
    # Input/output arguments
    parser.add_argument(
        "input",
        metavar="path",
        type=file_path,
        nargs="+",
        help="Report files from kir-mapper",
    )
    parser = add_out_altable_args(parser)

    parser.set_defaults(func=call_function)

    return parser


def call_function(args):
    """
    Main function to execute the kir-mapper report to allele table conversion.

    Args:
        args: Parsed command line arguments
    """
    reports = read_reports(args.input)
    reports = exclude_high_missings(reports, threshold=5)

    reports = parse_alleles(reports)

    all_consensus = dict()
    for sample, row in reports.iterrows():
        report = row.to_dict()

        parser = AlleleParser(
            gene_family=args.gene_family, config_file=args.config_file
        )
        gene = ConsensusGene(
            name=report["Gene"], calls=report["Calls"], allele_parser=parser
        )
        gene.set_consensus_settings(
            normalize_weight=False, max_support=len(report["Calls"])
        )

        consensus = gene.consensus_dict(min_support=0.6)
        consensus["original_calls"] = report["Calls"]

        all_consensus[sample] = consensus

    allele_table = pd.DataFrame.from_dict(all_consensus, orient="index")
    allele_table.index.name = "SampleID"

    gene = reports["Gene"].unique()[0]
    allele_table[[gene + "_1", gene + "_2"]
                 ] = allele_table["alleles"].apply(lambda x: pd.Series(x))

    alt = AlleleTable()
    alt.alleles = allele_table.drop(
        columns=["alleles", "original_calls", "coverage", "gene", "support"])
    alt.load_phenotype(args.phenotype)
    print(alt.alleles.head())

    if args.remove_pheno_zero:
        alt.remove_phenotype_zero()

    alt.to_csv(args.output)


def read_reports(files: List[str]) -> pd.DataFrame:
    all_dfs = list()
    for file in files:
        df = pd.read_csv(file, sep="\t", header=0, index_col=0)
        df["Gene"] = file.split(".")[0]
        all_dfs.append(df)

    return pd.concat(all_dfs)


def split_alleles(x: str):
    genotypes = x.split(";")

    # Add algorithm name (numbers) and filter alleles
    ret = dict()
    for idx, key in enumerate(genotypes):
        alleles = key.split("+")
        # Rename null to 000 alleles and unresolved to empty string
        alleles = [allele.replace("null", "000") for allele in alleles]
        alleles = [
            allele if "unresolved" not in allele else "" for allele in alleles]
        ret["kir-mapper" + str(idx)] = alleles
    return ret


def parse_alleles(df: pd.DataFrame) -> pd.DataFrame:
    alleles = df["Calls"].apply(split_alleles)
    df["Calls"] = alleles
    return df


def exclude_high_missings(df: pd.DataFrame, threshold: float = 5):
    # Get minimum missings
    df["MinMiss"] = df["Missings"].apply(get_min_number)

    # Filter alleles with high missings
    df = df[(df["MinMiss"] < threshold) | (df["MinMiss"].isna())]

    df.drop("MinMiss", axis=1)

    return df


def get_min_number(input: str):
    if not input or not isinstance(input, str):
        return input

    nums = [int(i) for i in input.split(";")]

    return min(nums)
