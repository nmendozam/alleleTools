"""
Interpret module for alleleTools.

This module provides tools for interpreting HLA allele data. So far, it includes
tools for visualizing individual characteristics of an HLA allele, helping
researchers understand the functional and clinical significance of allele
variants.
"""

from . import graph_pathogens, graph_phewas, plot_ikmb_coverage


def setup_parser(subparsers):
    """
    Set up the argument parser for the interpret subcommand.
    
    Creates a subparser for the 'interpret' command and registers all available
    interpretation tools as sub-subcommands.
    
    Args:
        subparsers: The subparsers object from the main argument parser.
        
    Returns:
        argparse.ArgumentParser: The configured parser for the interpret command.
    """
    parser = subparsers.add_parser(
        "plot",
        help="Get useful information about HLA alleles.",
        description="This commands provide some useful tools "
                    "for data interpretation. e.g. pleiotropy reports",
    )
    interpret_parser = parser.add_subparsers(
        dest="interpret_type",
        help="Available commands:",
        required=True,
    )

    graph_phewas.setup_parser(interpret_parser)
    graph_pathogens.setup_parser(interpret_parser)
    plot_ikmb_coverage.setup_parser(interpret_parser)

    return parser
