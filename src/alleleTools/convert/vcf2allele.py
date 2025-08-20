"""
VCF to Allele Table Conversion Module.

This module converts VCF (Variant Call Format) files
containing HLA or KIR allele data into allele tables. It supports various
output formats including pyHLA and PyPop compatible formats.

Usage:
    To generate the input file from imputation, run:
    # Extract only relevant alleles
    bcftools view --include 'ID~"HLA"' IMPUTED.vcf > HLA.vcf
    # Convert the extracted alleles to a table
    altools convert vcf2allele HLA.vcf --output out.alt

Author: Nicolás Mendoza Mejía (2023)
"""

from collections import defaultdict
from enum import Enum
from typing import List

import pandas as pd

from ..argtypes import file_path, output_path
from ..convert.vcf import VCF


def setup_parser(subparsers):
    """
    Set up the argument parser for the vcf2allele command.

    Args:
        subparsers: The subparsers object to add this command to.

    Returns:
        argparse.ArgumentParser: The configured parser for vcf2allele.
    """
    parser = subparsers.add_parser(
        name="vcf2allele",
        description="Convert vcf file to allele table",
        epilog="Author: Nicolás Mendoza Mejía (2023)",
    )
    # Input/output arguments
    parser.add_argument(
        "input",
        type=file_path,
        help="Input vcf file name",
    )
    parser.add_argument(
        "--phe",
        type=file_path,
        help="input phe file name (to add phenotype column)",
        default="",
    )
    parser.add_argument(
        "--output",
        type=output_path,
        help="name of the output file",
        default="output.alt",
    )
    # Allele format arguments
    parser.add_argument(
        "--rm-prefix",
        type=str,
        help="""
        removes prefix from allele names, usually its the gene name
        (KIR or HLA_)
        """,
        default="HLA_",
    )
    parser.add_argument(
        "--separator",
        type=str,
        help="separator to split gene name from allele name",
        default="*",
    )
    parser.add_argument(
        "--extensive_search",
        type=bool,
        help="""
        when no allele is imputed, look for the next most likely alleles
        """,
        default=False,
    )
    # Additional arguments
    parser.add_argument(
        "--output_header",
        action="store_true",
        help="output header with the gene names",
    )
    parser.add_argument(
        "--population",
        type=str,
        help="""
        If this is set, a colum with the population will be added at the
        beginning. This makes the output compatible with pyPop.
        """,
        default="",
    )

    parser.set_defaults(func=call_function)

    return parser


def call_function(args):
    """
    Main function to execute VCF to allele table conversion.

    This function orchestrates the conversion process by:
    1. Loading and preprocessing the VCF file
    2. Extracting genotype information
    3. Converting to allele format
    4. Adding phenotype data if provided
    5. Writing the output file

    Args:
        args: Parsed command line arguments containing:
            - input: Path to input VCF file
            - output: Path to output file
            - rm_prefix: Prefix to remove from allele names
            - separator: Separator between gene and allele names
            - extensive_search: Whether to perform extensive allele search
            - phe: Optional phenotype file path
            - output_header: Whether to include header in output
            - population: Population identifier for PyPop compatibility
    """
    vcf = VCF(args.input)
    vcf.remove_id_prefix(args.rm_prefix)
    format = vcf.get_format()

    # Get table with only samples as columns
    genotypes = vcf.samples_dataframe()

    true_alleles = _get_true_alleles(
        genotypes, format, args.extensive_search, args.separator
    )
    # sort the columns
    true_alleles = true_alleles.reindex(sorted(true_alleles.columns), axis=1)

    # add phenotype column if provided
    if args.phe:
        print("Warn: add phenotype (--phe) feature isn't tested")
        phe = pd.read_csv(args.phe, sep=" ", comment="##")
        phe.set_index("IID", inplace=True)
        true_alleles = true_alleles.join(phe.phenotype)
    else:
        true_alleles["phenotype"] = 0

    # move phenotype to the first column
    cols = true_alleles.columns.tolist()
    cols.insert(0, cols.pop(len(cols) - 1))

    true_alleles.reset_index(inplace=True)  # move sample id to column

    if args.population:
        true_alleles.insert(
            0, "population", args.population
        )  # add population column at the beginning

    true_alleles.to_csv(
        args.output,
        sep="\t",
        index=False,
        na_rep="NA",
        header=args.output_header
    )


class VCFalleles:
    """
    Class for processing and analyzing VCF allele data.

    This class handles the parsing and analysis of allele information from
    VCF format data, including genotype determination, resolution analysis,
    and ploidy-based filtering.

    Attributes:
        df (pd.DataFrame): Processed allele data with genotype information
                          indexed by gene and allele names.

    Args:
        alleles (pd.DataFrame): Raw allele data from VCF
        formats (pd.DataFrame): Format information from VCF header
    """

    def __init__(self, alleles: pd.DataFrame, formats: pd.DataFrame) -> None:
        self.df = self._parse_vcf_info(alleles, formats)

        self.df["is_homozygous"] = self.df.GT.str.contains(r"1\|1")

        genes = self._get_genes(self.df.index)
        self.df = self.df.set_index([genes, self.df.index]).rename_axis(
            ["gene", "allele"]
        )

    def _parse_vcf_info(self, alleles, formats):
        """
        Parse VCF info fields into structured data.

        Extracts genotype, dosage, and ploidy likelihoods from VCF format
        fields following the pattern: {GT}:{DS}:{AA},{AB},{BB}

        Args:
            alleles (pd.Series): Series with info fields from VCF file
            formats (list): List of format field names

        Returns:
            pd.DataFrame: Parsed data with GT as string and DS, AA, AB, BB as floats
        """
        # Extract the genotype, dosage and ploidy likelihoods
        info_df = alleles.str.split(":", expand=True)
        info_df.columns = formats

        float_cols = ["DS", "AA", "AB", "BB"]

        # Convert columns to the correct data type
        for col in set(formats).intersection(float_cols):
            info_df[col] = info_df[col].astype(float)

        return info_df

    def _find_high_res(self, allele_list: List[str]):
        """
        Identify high-resolution alleles from a list.

        Determines which alleles represent the highest resolution by checking
        if any allele is a substring of another (indicating lower resolution).
        For example, in ['A*02', 'A*02:01:01:01'], only 'A*02:01:01:01' is
        considered high resolution.

        Args:
            allele_list (List[str]): List of allele names to analyze

        Returns:
            List[bool]: Boolean list indicating which alleles are high resolution
        """
        high_res = list()
        for i, x in enumerate(allele_list):
            # cross comparison with every other allele
            comp = [x in y for y in allele_list]
            comp[i] = False  # remove self comparison
            high_res.append(not any(comp))

        return high_res

    def _get_genes(self, alleles: pd.Series):
        """
        Extract gene names from allele identifiers.

        Args:
            alleles (pd.Series): Series of allele identifiers

        Returns:
            pd.Series: Extracted gene names
        """
        return alleles.str.extract(r"([A-Z0-9]+)")[0]

    class _ploidy_filter(Enum):
        """
        Enumeration for filtering alleles by ploidy status.

        Each enum value contains a tuple of (mask_boolean, regex_pattern):
        - HOMOZYGOUS: Filters for homozygous genotypes (1|1)
        - HETEROZYGOUS: Filters for heterozygous genotypes (not 0|0)
        - UNCERTAIN: Filters for uncertain genotypes (0|0)
        """
        HOMOZYGOUS = (True, r"1\|1")
        HETEROZYGOUS = (False, r"0\|0")
        UNCERTAIN = (True, r"0\|0")

    def _get_ploidy_alleles(
        self, alleles: pd.DataFrame, filter: _ploidy_filter, n: int = 1
    ):
        """
        Get the highest resolution alleles that match a ploidy filter.

        Args:
            alleles (pd.DataFrame): DataFrame with allele codes as index and
                                  columns GT, DS, AB
            filter (_ploidy_filter): Filter criteria for ploidy selection
            n (int): Number of alleles to return (default: 1)

        Returns:
            pd.DataFrame: Filtered alleles sorted by score (DS + AB)
        """
        # 1. Apply filter
        mask, pattern = filter.value
        has_passed_filter = (
            alleles.GT.str.contains(pattern)
            if mask
            else ~alleles.GT.str.contains(pattern)
        )
        filtered = alleles.loc[has_passed_filter]
        # 2. Get high resolution alleles
        is_high_res = self._find_high_res(filtered.index.values.tolist())
        high_res = filtered.loc[is_high_res].copy()
        # 3. Get the n alleles with the highest dosages + AB probability
        if "AB" in high_res and "DS" in high_res:
            high_res["score"] = high_res["DS"] + high_res["AB"]
        else:
            high_res["score"] = 1

        return high_res.nlargest(1, columns="score")

    def _fill_ploidy_second_most_probable(self, heterozygous, alleles):
        """
        Fill ploidy with second most probable alleles when needed.

        When there aren't enough high-confidence alleles to fill the expected
        ploidy, this method finds the next most probable alleles.

        Args:
            heterozygous (pd.DataFrame): Current heterozygous alleles
            alleles (pd.DataFrame): All available alleles

        Returns:
            List[str]: Combined list of allele names including scores
        """
        # If there is not enough alleles to fill the ploidy
        # with the most probable alleles
        hetero_low = self._get_ploidy_alleles(
            alleles, self._ploidy_filter.UNCERTAIN, n=2 - len(heterozygous)
        )

        if not ("DS" in hetero_low and "AB" in hetero_low):
            return []

        hetero_low_str = hetero_low.apply(
            lambda row: f"{row.name}({row['DS']}/{row['AB']})", axis=1
        )

        return heterozygous.index.tolist() + [hetero_low_str.iloc[0]]

    def sort_and_fill(self, extensive=False):
        """
        Sort alleles by gene and determine final genotype calls.

        For each gene, determines whether the genotype is homozygous or
        heterozygous and selects the appropriate alleles based on confidence
        scores and resolution.

        Args:
            extensive (bool): Whether to perform extensive search for
                            low-confidence alleles (default: False)

        Returns:
            List[str]: List of selected allele names for all genes
        """
        results = list()
        genes = self.df.index.get_level_values("gene").unique()
        for gene in sorted(genes):
            to_append = list()
            alleles = self.df.loc[[gene]].reset_index(level=0)

            if any(alleles.is_homozygous):
                # If the gene is homozygous, get the allele
                homozygous = self._get_ploidy_alleles(
                    alleles, self._ploidy_filter.HOMOZYGOUS
                )
                to_append = homozygous.index.tolist() * 2
            else:  # is heterozygous
                # Get the two highest resolution alleles
                heterozygous = self._get_ploidy_alleles(
                    alleles, self._ploidy_filter.HETEROZYGOUS, n=2
                )
                if len(heterozygous == 2):
                    to_append = heterozygous.index.tolist()
                elif len(heterozygous) <= 1 and extensive:
                    to_append = self._fill_ploidy_second_most_probable(
                        heterozygous, alleles
                    )
                else:
                    to_append = []

            results.extend(to_append)

        return results


def _get_true_alleles(genotypes,
                      format,
                      extensive=False,
                      allele_separator="*"
                      ):
    """
    Extract true alleles from VCF genotype data for all samples.

    This function processes genotype data for all samples in the VCF file,
    determining the most likely allele calls for each gene and sample.

    Args:
        genotypes (pd.DataFrame): DataFrame with samples as columns and
                                allele IDs as index, containing genotype info
        format (list): List of format field names from VCF header
        extensive (bool): Whether to perform extensive search for uncertain
                        alleles (default: False)
        allele_separator (str): Character separating gene name from allele
                              name (default: "*")

    Returns:
        pd.DataFrame: DataFrame with samples as index and genes as columns,
                     containing the called alleles

    Raises:
        Warning: Prints warning if alleles are ignored due to separator mismatch
    """
    df = pd.DataFrame()
    ignored_samples = defaultdict(lambda: 0)
    last_ignored_allele = ""
    for sample in genotypes.columns:
        # Get the list of alleles for that column(sample)
        alleles = genotypes.loc[
            ~genotypes[sample].str.contains(r"0\|0:0:1,0,0", na=False), sample
        ]

        allele_list = VCFalleles(alleles, format).sort_and_fill(extensive)

        # Separate allele from gene name
        genes, alleles = (list(), list())
        for allele in allele_list:
            if allele_separator not in allele:
                ignored_samples[sample] += 1
                last_ignored_allele = allele
            else:
                gene, allele = allele.split(allele_separator)
                genes.append(gene)
                alleles.append(allele)

        if len(alleles) == 0:
            continue

        # Add _1 to duplicate genes
        genes_columns = list()
        for gene in genes:
            if gene in genes_columns:
                gene += "_1"
            genes_columns.append(gene)

        # Add the alleles to the data frame
        row = pd.DataFrame(
            allele_list,
            index=genes_columns,
            columns=[sample],
        ).transpose()
        df = pd.concat([df, row])

    if len(ignored_samples) > 0:
        print(
            f"""
            Warning: A total of {len(ignored_samples)} samples had ignored
            alleles. This may be due to a separator not matching the allele
            names, try changing the --separator argument. For example, the last
            ignored allele was '{last_ignored_allele}', but the expected
            --separator was '{allele_separator}'
            """
        )
    return df
