import json

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

    parser.set_defaults(func=call_function)

    return parser


def call_function(args):
    df = read_reports_asdf(args.input)

    ax = df.boxplot(column="coverage", by="gene")
    ax.set_ylim([0, 500])
    plt.show()


def read_reports_asdf(files: list) -> pd.DataFrame:
    """
    Parses a list of .json files as reports and returns
    a dataframe with all the genotyped genes per sample
    """
    reports = list()
    for file in files:
        j = read_json(file)
        report = Report(j)
        reports.extend(report.aslist())

    return pd.DataFrame(reports)
