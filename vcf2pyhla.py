import argparse
import re
from collections import defaultdict
from typing import List

import pandas as pd

"""
To generate the input file from the imputation run this command
> grep "#CHR" IMPUTED.vcf > HLA.vcf
> grep "HLA_" IMPUTED.vcf >> HLA.vcf
"""


def read_vcf(file_name):
    """
    This takes a vcf file data frame and returns
    a table were the row indexes are the allele names
    and the column names are the sample names
    """
    df = pd.read_csv(file_name, sep="\t", on_bad_lines="warn")
    # Use alleles as index
    df["ID"] = df["ID"].str.replace("HLA_", "")
    df.set_index("ID", inplace=True)
    # Drop non sample columns
    df.drop(
        ["#CHROM", "POS", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"],
        axis=1,
        inplace=True,
    )
    return df


class AlleleList:
    def __init__(self, alleles: pd.DataFrame) -> None:
        self.df = self._parse_vcf_info(alleles)

        self.df["is_homozygous"] = self.df.GT.str.contains("1\|1")

        genes = self._get_genes(self.df.index)
        self.df = self.df.set_index([genes, self.df.index]).rename_axis(
            ["gene", "allele"]
        )

    def _parse_vcf_info(self, alleles):
        """
        This function takes a pandas series with info fields from a vcf file
        and the allele codes as index. It returns a data frame with the
        genotype, dosage and ploidy likelihoods. The pattern is:

        {GT}:{DS}:{AA},{AB},{BB}

        GT is a string and the other four are floats
        """
        # Extract the genotype, dosage and ploidy likelihoods
        pattern = r"(?P<GT>[^:]*):(?P<DS>[^:]*):(?P<AA>[^,]*),(?P<AB>[^,]*),(?P<BB>.*)"
        info_df = alleles.str.extract(pattern)

        # Convert columns to the correct data type
        for col in ["DS", "AA", "AB", "BB"]:
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
        return alleles.str.extract(r"([A-Z0-1]+?)\*")[0]

    def sort_and_fill(self):
        results = list()
        genes = self.df.index.get_level_values("gene").unique()
        for gene in sorted(genes):
            to_append = list()
            alleles = self.df.loc[[gene]].reset_index(level=0)

            if any(alleles.is_homozygous):
                # 1. get only homozygous flag 1|1
                homozygous = alleles.loc[alleles.is_homozygous]
                # 2. get the highest resolution
                homozygous = homozygous.loc[
                    self._find_high_res(homozygous.index.values.tolist())
                ]
                to_append = homozygous.index.tolist() * 2
            else:  # is hetereozygous
                # 1. remove all 0|0
                heterozygous = alleles.loc[~alleles.GT.str.contains("0\|0")]
                # 2. Get high resolution alleles
                heterozygous = heterozygous.loc[
                    self._find_high_res(heterozygous.index.values.tolist())
                ]
                # 3. Get the two with the highest dosages + AB probability
                heterozygous["DS_AB"] = heterozygous["DS"] + heterozygous["AB"]
                heterozygous = heterozygous.nlargest(2, columns="DS_AB")
                if len(heterozygous) == 1:
                    # 4. If there is only one high resolution allele, find the next highest resolution 0|0
                    hetero_low = alleles.loc[alleles.GT.str.contains("0\|0")]
                    hetero_low = hetero_low.loc[
                        self._find_high_res(hetero_low.index.values.tolist())
                    ]
                    # 5. Get the highest dosage + AB probability from posible alleles
                    hetero_low["DS_AB"] = hetero_low["DS"] + hetero_low["AB"]
                    hetero_low = hetero_low.nlargest(1, columns="DS_AB")

                    # from str to output dosage and AB
                    hetero_low_str = hetero_low.apply(
                        lambda row: f"{row.name}({row['DS']}/{row['AB']})", axis=1
                    )

                    to_append = heterozygous.index.tolist() + [hetero_low_str.iloc[0]]
                elif len(heterozygous == 2):
                    to_append = heterozygous.index.tolist()
                else:
                    to_append = ["NA", "NA"]

            results.extend(to_append)

        return results


def get_true_alleles(vcf):
    vcf = read_vcf(vcf)
    true_alleles = dict()
    for sample in vcf.columns:
        # Get the list of alleles for that column(sample)
        alleles = vcf.loc[~vcf[sample].str.contains("0\|0:0:1,0,0", na=False), sample]

        true_alleles[sample] = AlleleList(alleles).sort_and_fill()
    return true_alleles


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convert vcf file to pyhla and output to std out"
    )
    parser.add_argument("vcf", type=str, help="input vcf file name")
    parser.add_argument("--phe", type=str, help="input phe file name", default="")

    args = parser.parse_args()

    if args.phe:
        phe = pd.read_csv(args.phe, sep=" ", comment="##")
        phe.set_index("IID", inplace=True)

    true_alleles = get_true_alleles(args.vcf)
    for key in true_alleles:
        phenotype = 1
        if args.phe:
            phenotype = int(phe.loc[[key]].LLI)
        print(key, phenotype, "\t".join(true_alleles[key]), sep="\t")
