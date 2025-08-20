"""
Convert module for alleleTools.

This module converts allele data between different
file formats including VCF, CSV, and custom allele table formats. It supports
conversions from various sources including UK Biobank data and standard
genomic formats.
"""

from . import allele2vcf, ukb2allele, vcf2allele


def setup_parser(subparsers):
    """
    Set up the argument parser for the convert subcommand.

    Creates a subparser for the 'convert' command and registers all available
    conversion types as sub-subcommands.

    Args:
        subparsers: The subparsers object from the main argument parser.

    Returns:
        argparse.ArgumentParser: The configured parser for the convert command.
    """
    parser = subparsers.add_parser(
        "convert",
        help="Convert allele data between different formats.",
        description="""
        This command facilitates the conversion of allele data from allele
        tables to other formats (e.g., VCF or CSV). It supports various
        input/output formats and allows for optional validation during
        conversion. The validations are specific to each conversion type check
        the --help of each conversion type for further details.
        """,
    )
    convert_parser = parser.add_subparsers(
        dest="convert_type",
        help="Specify the type of conversion to perform.",
        required=True,
    )

    # Register toolconvertion tools
    vcf2allele.setup_parser(convert_parser)
    allele2vcf.setup_parser(convert_parser)
    ukb2allele.setup_parser(convert_parser)

    return parser
