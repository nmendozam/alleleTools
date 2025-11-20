"""
Convert module for alleleTools.

This module converts allele data between different
file formats including VCF, CSV, and custom allele table formats. It supports
conversions from various sources including UK Biobank data and standard
genomic formats.
"""

from alleleTools.format import allele_resolution

from . import (allele2vcf, from_ikmb_hla, hla_group, immuannot_report,
               kir_mapper, ukb2allele, vcf2allele)


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
        "format",
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

    # Register conversion tools
    vcf2allele.setup_parser(convert_parser)
    allele2vcf.setup_parser(convert_parser)
    ukb2allele.setup_parser(convert_parser)
    kir_mapper.setup_parser(convert_parser)
    from_ikmb_hla.setup_parser(convert_parser)
    allele_resolution.setup_parser(convert_parser)
    immuannot_report.setup_parser(convert_parser)
    hla_group.setup_parser(convert_parser)

    return parser
