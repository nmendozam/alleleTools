"""
Genotype module for alleleTools.

This module has some tools useful during the genotyping process of polymorphic
genes. It includes consensus algorithms for combining results from multiple
genotyping methods and quality assessment tools.
"""

from . import consensus


def setup_parser(subparsers):
    """
    Set up the argument parser for the genotype subcommand.
    
    Creates a subparser for the 'genotype' command and registers all available
    genotyping tools as sub-subcommands.
    
    Args:
        subparsers: The subparsers object from the main argument parser.
        
    Returns:
        argparse.ArgumentParser: The configured parser for the genotype command.
    """
    parser = subparsers.add_parser(
        "genotype",
        help="Useful commands to use when genotyping polymorphic genes.",
        description="This command has some useful algorithms to use when genotyping polymorphic genes.",
    )
    genotype_parser = parser.add_subparsers(
        dest="genotype_type",
        help="Specify the operation to perform.",
        required=True,
    )

    consensus.setup_parser(genotype_parser)

    return parser
