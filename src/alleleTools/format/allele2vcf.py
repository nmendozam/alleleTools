"""
Allele Table to VCF Conversion Module.

This module provides functionality to convert allele table data into VCF
(Variant Call Format) files. It handles the transformation of genotype data
from tabular format to standard VCF format, including proper diploid notation
and genomic coordinate mapping.

Author: Nicolás Mendoza Mejía (2023)
"""

import pandas as pd

from ..argtypes import file_path
from ..utils.assets import get_asset_path


def setup_parser(subparsers):
    """
    Set up the argument parser for the allele2vcf command.

    Args:
        subparsers: The subparsers object to add this command to.

    Returns:
        argparse.ArgumentParser: The configured parser for allele2vcf.
    """
    parser = subparsers.add_parser(
        name="to_vcf",
        help="Convert allele table to vcf",
        description="Convert allele table to vcf",
        epilog="Author: Nicolás Mendoza Mejía (2023)",
    )
    parser.add_argument(
        "input",
        type=str,
        help="path to the input file with allele data in allele table format",
    )
    parser.add_argument(
        "--loci_file",
        type=file_path,
        help="""
        path to the file with gene loci information, alternatively you could
        specify --gene_cluster
        """,
    )
    parser.add_argument(
        "--gene_cluster",
        type=str,
        choices=["HLA", "KIR"],
        help="""
        name of the gene cluster (hla or kir), alternatively you could provide
        a --loci_file
        """,
    )
    parser.add_argument(
        "--vcf",
        type=file_path,
        help="VCF file to append the alleles to",
        default="file.vcf",
        required=True,
    )
    parser.add_argument(
        "--field_separator",
        type=str,
        help="character separating the fields inside input default is tab",
        default="\t",
    )

    parser.set_defaults(func=call_function)

    return parser


def call_function(args):
    """
    Main function to execute allele table to VCF conversion.

    This function orchestrates the conversion process by:
    1. Validating input parameters
    2. Loading and processing allele table data
    3. Converting to VCF format with proper diploid notation
    4. Mapping alleles to genomic coordinates
    5. Appending results to the target VCF file

    Args:
        args: Parsed command line arguments containing:
            - input: Path to input allele table file
            - loci_file: Path to gene loci information file
            - gene_cluster: Alternative gene cluster specification
            - vcf: Path to target VCF file for appending
            - field_separator: Field separator for input file (default: tab)

    Raises:
        SystemExit: If neither gene_cluster nor loci_file is provided
    """
    if not (args.gene_cluster or args.loci_file):
        print(
            "Error: either --gene_cluster or --loci_file must be provided."
            "Use -h or --help to see more details."
        )
        exit(1)
    
    if not args.loci_file:
        if args.gene_cluster == "HLA":
            args.loci_file = get_asset_path("gene_table.tsv")
        elif args.gene_cluster == "KIR":
            args.loci_file = get_asset_path("gene_table_kir.tsv")
        else:
            print(
                "Error: --gene_cluster must be either 'HLA' or 'KIR'."
                "Use -h or --help to see more details."
            )
            exit(1)

    genotypes = pd.read_csv(args.input, sep=args.field_separator)
    gene_loci = pd.read_csv(args.loci_file, sep="\t")

    pairs = _gene_pairs(genotypes.columns)

    # By this point we have the vcf file without the leading columns.
    # First we pivot the table to have the samples as columns and the alleles
    # as rows. Then we merge both gene alleles into a single row by applying
    # the diploid_notation function. e.g.:
    # ID  SAMPLE_ID ...
    pre_vcf_alleles = pd.DataFrame()
    for pair in pairs:
        # Get the presence/absence of each allele in the samples.
        pairA, pairB = pair
        pivotA = genotypes.pivot_table(
            index=pairA, columns="id", values=pairB, aggfunc="sum"
        ).notna()
        pivotB = genotypes.pivot_table(
            index=pairB, columns="id", values=pairA, aggfunc="sum"
        ).notna()
        allele_codes = _diploid_notation(pivotA, pivotB)
        # add column with gene name
        allele_codes["gene"] = pairA
        pre_vcf_alleles = pd.concat([pre_vcf_alleles, allele_codes])

    # Now we add the leading columns and sort the samples (also in columns)
    # to match the base vcf file order.
    # e.g.:
    # CHROM  POS ID  REF  ALT  QUAL  FILTER  INFO  FORMAT  SAMPLE_ID. ...
    # column start has chr:pos, divide that into two columns.
    gene_loci[["CHROM", "POS"]] = gene_loci["start"].str.split(
        ":", expand=True)
    # cross gene from pre_vcf_alleles with gene_loci to get the position of
    # each allele.
    vcf_alleles = pre_vcf_alleles.merge(gene_loci, how="left", on="gene")
    vcf_alleles = vcf_alleles[vcf_alleles["POS"].notna()]
    vcf_alleles = vcf_alleles.assign(
        REF="A", ALT="T", QUAL=".", FILTER="PASS", INFO=".", FORMAT="GT"
    )
    # Remove gene columns and rename
    vcf_alleles = vcf_alleles.drop(columns=["gene", "start"])

    # Sort columns to match the base vcf file order.
    vcf_col = _get_vcf_columns(args.vcf)
    vcf_col = [x for x in vcf_col if x in vcf_alleles.columns]
    vcf_alleles = vcf_alleles[vcf_col]
    vcf_alleles.fillna("0|0", inplace=True)

    if len(vcf_alleles) == 0:
        print("WARNING: No alleles are being added to the VCF file.")
        print("Check that the gene names in the input file (column names) match those in the loci file.")
    else:
        print(f"Appending {len(vcf_alleles)} alleles to {args.vcf}")
    with open(args.vcf, "a") as f:
        f.write(vcf_alleles.to_csv(index=False, sep="\t", header=False))


def _diploid_notation(pairA, pairB):
    """
    Convert allele presence/absence data to VCF diploid notation.

    Takes two boolean DataFrames representing allele presence for each gene
    copy and converts them to standard VCF genotype notation (0|0, 0|1, 1|0,
    1|1).

    Args:
        pairA (pd.DataFrame): Boolean DataFrame for first allele copy
        pairB (pd.DataFrame): Boolean DataFrame for second allele copy

    Returns:
        pd.DataFrame: DataFrame with allele IDs and corresponding VCF
            genotype codes

    Example:
        pairA=True, pairB=False -> "1|0"
        pairA=True, pairB=True -> "1|1"
    """
    # assign a code to each allele combination.
    pairA = pairA.astype(int).rename_axis("id")
    pairB = pairB.astype(int).mul(2).rename_axis("id")
    group = pd.concat([pairA, pairB]).groupby("id").sum()
    group = group.replace({0: "0|0", 1: "1|0", 2: "0|1", 3: "1|1"})
    # Add gene name to the index
    group.index = group.index.str.replace("*", "_")
    group.index.name = "ID"
    return group.reset_index()


def _gene_pairs(lst):
    """
    Extract gene pairs from column names following the gene/gene.1 convention.

    Identifies paired gene columns where one column represents the first allele
    (e.g., "HLA-A") and another represents the second allele (e.g., "HLA-A.1").

    Args:
        lst (list): List of column names from the allele table

    Returns:
        list: List of tuples containing paired gene column names
              [(gene, gene.1), ...]

    Example:
        >>> _gene_pairs(["HLA-A", "HLA-A.1", "HLA-B", "HLA-B.1"])
        [["HLA-A", "HLA-A.1"], ["HLA-B", "HLA-B.1"]]
    """
    # Get elements in the list with *.1
    gene_1 = [x for x in lst if x.endswith(".1")]
    # Get possible elements without *.1
    posible = [x.replace(".1", "") for x in gene_1]
    # Get elements in the list without *.1 that have a pair
    gene = [x for x in lst if x not in gene_1 and x in posible]
    # Remove elements without a pair in gene_1
    gene_1 = [x for x in gene_1 if x.replace(".1", "") in gene]

    return [[g, g1] for g, g1 in zip(gene, gene_1)]


def _get_vcf_columns(vcf_file):
    """
    Extract column names from a VCF file header.

    Reads the VCF file to find the #CHROM header line and extracts
    the column names for proper column ordering in the output.

    Args:
        vcf_file (str): Path to the VCF file

    Returns:
        list: List of column names from the VCF header

    Example:
        >>> _get_vcf_columns("sample.vcf")
        ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
        "FORMAT", "SAMPLE1", "SAMPLE2"]
    """
    # Read only the line that starts with #CHROM
    with open(vcf_file, "r") as f:
        line = f.readline()
        while not line.startswith("#CHROM"):
            line = f.readline()
    # Remove leading # and \n, then split by tab.
    return line[1:].strip().split("\t")  # [9:]
    return line[1:].strip().split("\t")  # [9:]
    return line[1:].strip().split("\t")  # [9:]
