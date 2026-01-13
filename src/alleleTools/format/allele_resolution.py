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

from .alleleTable import AlleleTable


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
    parser.add_argument(
        "--max_miss",
        type=int,
        help="maximum missing alleles allowed per sample",
        default=0,
    )

    parser.set_defaults(func=call_function)


def call_function(args):
    """
    Normalize allele resolutions with the provided arguments.
    """
    parser = AlleleParser(gene_family=args.gene_family, config_file=args.config_file)

    alt = AlleleParsedTable.open(args.input)
    alt.parse_alleles(allele_parser=parser)
    alt.normalize_resolution(resolution=args.resolution)

    # filter out samples with too many missing alleles
    if args.max_miss > 0:
        alt = alt.convert_to_altable()
        mask = alt.alleles.isnull().sum(axis=1) < args.max_miss
        alt.alleles = alt.alleles[mask]
        print(
            f"Filtered out {(~mask).sum()} samples with more than {args.max_miss} missing alleles"
        )

    alt.to_csv(args.output)


class AlleleParsedTable(AlleleTable):
    @classmethod
    def open(cls, filename: str, sep: str = "\t") -> "AlleleParsedTable":
        base = AlleleTable.open(filename, sep)
        inst = cls()
        inst.alleles = base.alleles
        inst.phenotype = base.phenotype
        inst.covariates = base.covariates
        return inst

    def normalize_resolution(self, resolution: int) -> "AlleleParsedTable":
        self.alleles = self.alleles.map(lambda x: x.truncate(resolution))
        return self

    def parse_alleles(self, allele_parser) -> "AlleleParsedTable":
        df = self.alleles.copy()
        df.fillna("", inplace=True)
        df = df.map(lambda x: allele_parser.parse(x))
        self.alleles = df
        return self

    def convert_to_altable(self) -> AlleleTable:
        alt = AlleleTable()
        alt.alleles = self._alleles_as_str_()
        alt.phenotype = self.phenotype.copy()
        alt.covariates = self.covariates.copy()
        return alt

    def to_csv(
            self, filename: str, header: bool = True, population: str = ""
    ):
        """
        Export the allele table to a CSV file.

        Args:
            filename (str): The name of the output CSV file.
            header (bool): Flag to store the file with column names or not
            population (str): Adds an extra column in the position left to
                phenotype with a population name. Currently, only one
                population per allele table is supported.
        """
        # Convert alleles back to string
        self.alleles = self.alleles.astype(str)
        super().to_csv(filename, header, population)
