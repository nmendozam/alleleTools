import argparse
from enum import Enum
from typing import List

import pandas as pd

"""
To generate the input file from the imputation run this command
> bcftools view --include 'ID~"HLA"' IMPUTED.vcf > HLA.vcf
"""


def read_vcf(file_name, prefix):
    """
    This takes a vcf file data frame and returns
    a table were the row indexes are the allele names
    and the column names are the sample names
    """
    last_pos = 0
    with open(file_name, "r") as f:
        # Skip the header
        while True:
            line = f.readline()
            if line.startswith("#CHROM"):
                break
            last_pos = f.tell()
        # Read the rest of the file
        f.seek(last_pos)
        df = pd.read_csv(f, sep="\t", on_bad_lines="warn")
        # Use alleles as index
        df["ID"] = df["ID"].str.replace(prefix, "")
        df.set_index("ID", inplace=True)

        # Get the format ofeeach allele
        format = df["FORMAT"].str.split(":", expand=True).iloc[0].tolist()

        # Drop non sample columns
        df.drop(
            ["#CHROM", "POS", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"],
            axis=1,
            inplace=True,
        )
        return df, format


class AlleleList:
    def __init__(self, alleles: pd.DataFrame, formats: pd.DataFrame) -> None:
        self.df = self._parse_vcf_info(alleles, formats)

        self.df["is_homozygous"] = self.df.GT.str.contains("1\|1")

        genes = self._get_genes(self.df.index)
        self.df = self.df.set_index([genes, self.df.index]).rename_axis(
            ["gene", "allele"]
        )

    def _parse_vcf_info(self, alleles, formats):
        """
        This function takes a pandas series with info fields from a vcf file
        and the allele codes as index. It returns a data frame with the
        genotype, dosage and ploidy likelihoods. The pattern is:

        {GT}:{DS}:{AA},{AB},{BB}

        GT is a string and the other four are floats
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
        Checks wether the alleles is covered by a higher resolution one
        or not. Because the imputation returns multiple levels of resolutions
        eg.  ['A*02', 'A*02:01:01:01' ] returns [False, True]
        """
        high_res = list()
        for i, x in enumerate(allele_list):
            # cross comparison with every other allele
            comp = [x in y for y in allele_list]
            comp[i] = False  # remove self comparison
            high_res.append(not any(comp))

        return high_res

    def _get_genes(self, alleles: pd.Series):
        return alleles.str.extract(r"([A-Z0-9]+)")[0]

    # enum of filter to get alleles by ploidy
    class _ploidy_filter(Enum):
        HOMOZYGOUS = (True, "1\|1")
        HETEROZYGOUS = (False, "0\|0")
        UNCERTAIN = (True, "0\|0")

    def _get_ploidy_alleles(
        self, alleles: pd.DataFrame, filter: _ploidy_filter, n: int = 1
    ):
        """
        Gets the n highest resolution alleles that have passed the filter.
        args:
            alleles: data frame with the allele codes as index and the
                columns GT, DS, AB
            filter: tuple of mask (get values that do or don't match the
                pattern) and pattern (substring to match)
            n: int number of alleles to return
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


def get_true_alleles(genotypes, format, extensive=False, allele_separator="*"):
    df = pd.DataFrame()
    for sample in genotypes.columns:
        # Get the list of alleles for that column(sample)
        alleles = genotypes.loc[
            ~genotypes[sample].str.contains("0\|0:0:1,0,0", na=False), sample
        ]

        allele_list = AlleleList(alleles, format).sort_and_fill(extensive)
        genes, alleles = zip(*[x.split(allele_separator) for x in allele_list])

        # Add _1 to duplicate genes
        genes_columns = list()
        for gene in genes:
            if gene in genes_columns:
                gene += "_1"
            genes_columns.append(gene)

        # Add the alleles to the data frame
        row = pd.DataFrame(
            allele_list, index=genes_columns, columns=[sample]
        ).transpose()
        df = pd.concat([df, row])

    return df


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convert vcf file to allele table and output to std out"
    )
    ## Input/output arguments
    parser.add_argument("vcf", type=str, help="input vcf file name")
    parser.add_argument(
        "--phe",
        type=str,
        help="input phe file name (to add phenotype column)",
        default="",
    )
    parser.add_argument(
        "--out", type=str, help="name of the output file", default="output.pyhla"
    )
    ## Allele format arguments
    parser.add_argument(
        "--prefix",
        type=str,
        help="prefix to remove from allele names(genename)",
        default="HLA_",
    )
    parser.add_argument(
        "--separator",
        type=str,
        help="separator to split gene name from allele name",
        default="*",
    )
    parser.add_argument(
        "--extensive",
        type=bool,
        help="when no allele is imputed, look for the next most likely alleles",
        default=False,
    )
    ## Additional arguments
    parser.add_argument(
        "-output_header",
        action="store_true",
        help="output header with the gene names",
    )
    parser.add_argument(
        "--population",
        type=str,
        help="""If this is set, a colum with the population will be added at the beginning.
                This makes the output compatible with pyPop""",
        default="",
    )

    args = parser.parse_args()

    genotypes, format = read_vcf(args.vcf, args.prefix)
    true_alleles = get_true_alleles(genotypes, format, args.extensive, args.separator)
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
        args.out, sep="\t", index=False, na_rep="NA", header=args.output_header
    )
