"""
Refactor module for alleleTools.

This module provides tools for refactoring and normalizing allele data,
including resolution standardization. These tools help ensure consistency across
different allele datasets and facilitate downstream analysis.
"""

from . import allele_resolution_wrapper


def setup_parser(subparsers):
    """
    Set up the argument parser for the refactor subcommand.
    
    Creates a subparser for the 'refactor' command and registers all available
    refactoring tools as sub-subcommands.
    
    Args:
        subparsers: The subparsers object from the main argument parser.
        
    Returns:
        argparse.ArgumentParser: The configured parser for the refactor command.
        
    Note:
        Some refactoring tools are currently implemented as shell scripts
        and are not yet integrated into the Python CLI interface.
    """
    parser = subparsers.add_parser(
        "refactor",
        help="Refactor allele table data.",
        description="This command facilitates usual refactoring of allele table data.",
    )
    refactor_parser = parser.add_subparsers(
        dest="refactor_type",
        help="Specify the type of refactoring to perform.",
        required=True,
    )

    allele_resolution_wrapper.setup_parser(refactor_parser)

    return parser
