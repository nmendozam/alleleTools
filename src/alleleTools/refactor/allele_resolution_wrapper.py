"""
Allele Resolution Normalization CLI Wrapper

This module provides a command-line interface for normalizing allele resolutions.
All alleles in the input file will be normalized to the specified resolution.
Resolutions are given in fields (e.g., 'one', 'two', 'three').

Example:
    input_file:
        HLA-A*01:01:01:01	HLA-B*01:01:01:01
        HLA-A*02:01:01:02	HLA-B*01:01:01:02
    command:
        altools allele_resolution two input_file output_file
    output_file:
        HLA-A*01:01	HLA-B*01:01
        HLA-A*02:01	HLA-B*01:01

Author: Nicolás Mendoza Mejía (2025)
"""
import os


def setup_parser(subparsers):
    parser = subparsers.add_parser(
        name="allele_resolution",
        help="Normalize allele resolutions",
        description="Normalize allele resolutions",
        epilog="Author: Nicolás Mendoza Mejía (2025)",
    )
    parser.add_argument(
        "resolution",
        help="The resolution to normalize (e.g., 'one', 'two', 'three')",
        type=str,
        choices=["one", "two", "three"],
        default="three",
    )
    parser.add_argument(
        "input",
        type=str,
        help="Path to the input file containing allele resolutions",
    )
    parser.add_argument(
        "output",
        type=str,
        help="Path to the output file where normalized allele resolutions will be saved",
    )
    parser.add_argument(
        "--prefix",
        type=str,
        help="Symbol separating gene name from allele number. For HLA-A*01:02 is * (default: '\\*')",
        default="\\*",
    )

    parser.set_defaults(func=call_function)


def call_function(args):
    """
    Normalize allele resolutions with the provided arguments.
    """

    # Get the base directory of the current script
    base_dir = os.path.dirname(os.path.abspath(__file__))
    cmd = ["bash", os.path.join(base_dir, "allele_resolution.sh")]

    # Add the arguments provided by the user
    cmd.append(args.resolution)
    cmd.append(args.input)
    cmd.append(args.output)
    cmd.append(args.prefix)