import json

from alleleTools.allele import AlleleParser
import matplotlib.pyplot as plt
import pandas as pd

from .ikmb_report import Report, read_json


def setup_parser(subparsers):
    """
    Set up the argument parser for the consensus command.

    Args:
        subparsers: The subparsers object to add this command to.

    Returns:
        argparse.ArgumentParser: The configured parser for consensus.
    """
    parser = subparsers.add_parser(
        name="plot_ikmb_coverage",
        description="""
        This program plots the coverage of HLA genes from IKMB reports.
        """,
        epilog="Author: Nicolás Mendoza Mejía (2023)",
    )
    parser.add_argument(
        "input",
        metavar="path",
        type=str,
        nargs="+",
        help="JSON files with HLA genotyping reports from the IKMB pipeline",
    )
    parser.add_argument(
        "--gene_family",
        type=str,
        help="Gene family to plot (e.g., 'hla')",
        default="hla",
    )
    parser.add_argument(
        "--config_file",
        type=str,
        help="Path to a custom allele parsing configuration file",
        default="",
    )

    parser.set_defaults(func=call_function)

    return parser


def call_function(args):
    parser = AlleleParser(gene_family=args.gene_family, config_file=args.config_file)
    df = read_reports_asdf(args.input, allele_parser=parser)

    ax = df.boxplot(column="coverage", by="gene")
    ax.set_ylim([0, 500])
    plt.show()


def read_reports_asdf(files: list, allele_parser: AlleleParser) -> pd.DataFrame:
    """
    Parses a list of .json files as reports and returns
    a dataframe with all the genotyped genes per sample
    """
    reports = list()
    for file in files:
        j = read_json(file)
        report = Report(j, allele_parser=allele_parser)
        reports.extend(report.aslist())

    return pd.DataFrame(reports)
