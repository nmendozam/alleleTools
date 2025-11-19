"""
Custom argument types for argparse validation. It validates file paths,
output paths and specific file formats.
"""

import argparse
import os
from pathlib import Path


def path(in_string: str) -> str:
    """
    Validate that a given path exists.

    Args:
        in_string (str): The path string to validate.

    Returns:
        str: The validated path string.

    Raises:
        argparse.ArgumentTypeError: If the path does not exist.

    Example:
        >>> path("/existing/path")
        "/existing/path"
    """
    if not os.path.exists(in_string):
        raise argparse.ArgumentTypeError(f"{in_string} does not exist.")
    return in_string


def output_path(in_string: str) -> str:
    """
    Validate that the output directory exists and the file doesn't already
    exist.

    Args:
        in_string (str): The output file path to validate.

    Returns:
        str: The validated output path string.

    Raises:
        argparse.ArgumentTypeError: If the parent directory doesn't exist
                                  or if the file already exists.

    Example:
        >>> output_path("/existing/dir/new_file.txt")
        "/existing/dir/new_file.txt"
    """
    basepath = Path(in_string).parent
    if not basepath.exists():
        raise argparse.ArgumentTypeError(
            f"the directory {basepath} does not exist.")
    elif os.path.isfile(in_string):
        raise argparse.ArgumentTypeError(f"{in_string} already exists!")
    return in_string


def file_path(in_string: str) -> str:
    """
    Validate that a file exists and is actually a file (not a directory).

    Args:
        in_string (str): The file path to validate.

    Returns:
        str: The validated file path string.

    Raises:
        argparse.ArgumentTypeError: If the path doesn't exist or is not a file.

    Example:
        >>> file_path("/path/to/existing_file.txt")
        "/path/to/existing_file.txt"
    """
    # Default parameters are empty, so we don't complain
    if not in_string:
        return in_string

    file = path(in_string)
    if not os.path.isfile(file):
        raise argparse.ArgumentTypeError(
            f"The file '{in_string}' is not a valid file")
    return file


def csv_file(in_string: str) -> str:
    """
    Validate that a CSV file exists and has the correct extension.

    Args:
        in_string (str): The CSV file path to validate.

    Returns:
        str: The validated CSV file path string.

    Raises:
        argparse.ArgumentTypeError: If the file doesn't exist, is not a file,
                                  or doesn't have a .csv extension.

    Example:
        >>> csv_file("/path/to/data.csv")
        "/path/to/data.csv"
    """
    file = file_path(in_string)
    filename, file_extension = os.path.splitext(file)
    if not file_extension.lower() == '.csv':
        raise argparse.ArgumentTypeError(
            f"{in_string} is not a valid CSV file")
    return file


def add_out_altable_args(parser):
    parser.add_argument(
        "--output",
        type=output_path,
        help="name of the output file",
        default="output.alt",
    )
    parser.add_argument(
        "--phenotype",
        type=str,
        help="""
        ssv file with 6 columns: eid, fid, ... , Sex, Pheno. No headers and
        space separated. The column Pheno (last column) will be included as
        phenotype in the output file.
        """,
        default="",
    )

    # Additional arguments
    parser.add_argument(
        "--remove_pheno_zero",
        action="store_true",
        help="Remove individuals with phenotype 0 from the output",
        default=False,
    )
    parser.add_argument(
        "--gene_family",
        type=str,
        help="Specify the gene family e.i. 'hla', 'kir'",
        default="kir",
    )
    parser.add_argument(
        "--config_file",
        type=file_path,
        help="Path to a custom allele parsing configuration file",
        default="",
    )
    return parser
