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

from alleleTools.argtypes import file_path
import pandas as pd

from alleleTools.allele import AlleleParser

from ..convert.alleleTable import AlleleTable


def setup_parser(subparsers):
    parser = subparsers.add_parser(
        name="allele_resolution",
        help="Normalize allele resolutions",
        description="Normalize allele resolutions",
        epilog="Author: Nicolás Mendoza Mejía (2025)",
    )
    parser.add_argument(
        "input",
        type=str,
        help="Path to the input file containing allele resolutions",
    )
    parser.add_argument(
        "--output",
        type=str,
        help="Path to the output file where normalized allele resolutions will be saved",
    )
    parser.add_argument(
        "--resolution",
        help="The resolution to normalize (e.g., 'one', 'two', 'three')",
        type=int,
        choices=[1, 2, 3],
        default=3,
    )
    parser.add_argument(
        "--gene_family",
        type=str,
        help="Specify the gene family e.i. 'hla', 'kir'",
        default="hla",
    )
    parser.add_argument(
        "--config_file",
        type=file_path,
        help="Path to a custom allele parsing configuration file",
        default="",
    )


    parser.set_defaults(func=call_function)


def call_function(args):
    """
    Normalize allele resolutions with the provided arguments.
    """
    parser = AlleleParser(gene_family=args.gene_family, config_file=args.config_file)
    alt = open_allele_table(args.input, allele_parser=parser)
    alt = normalize_resolution(alt, resolution=args.resolution)
    alt.to_csv(args.output)


def open_allele_table(input: str, allele_parser: AlleleParser, fields_separator: str = "\t") -> AlleleTable:
    alt = AlleleTable()
    df = pd.read_csv(input, sep=fields_separator)

    df.set_index("sample", inplace=True)

    alt.phenotype = df.pop("phenotype")
    alt.alleles = df

    alt = parse_allele_table(alt, allele_parser=allele_parser)

    return alt


def normalize_resolution(alt: AlleleTable, resolution: int) -> AlleleTable:
    alt.alleles = alt.alleles.map(lambda x: x.truncate(resolution))
    return alt


def parse_allele_table(alt: AlleleTable, allele_parser: AlleleParser) -> AlleleTable:
    df = alt.alleles.copy()
    df.fillna("", inplace=True)
    df = df.map(lambda x: allele_parser.parse(x))
    alt.alleles = df
    return alt
