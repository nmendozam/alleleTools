import argparse
from enum import Enum
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

    # enum of filter to get alleles by ploidy
    class _ploidy_filter(Enum):
        HOMOZYGOUS = (True, "1\|1")
        HETEROZYGOUS = (False, "0\|0")
        UNCERTAIN = (True, "0\|0")

    def _get_ploidy_alleles(self, alleles, filter, n=1):
        # 1. Apply filter
        mask, pattern = filter
        has_passed_filter = (
            alleles.GT.str.contains(pattern)
            if mask
            else ~alleles.GT.str.contains(pattern)
        )
        filtered = alleles.loc[has_passed_filter]
        # 2. Get high resolution alleles
        is_high_res = self._find_high_res(filtered.index.values.tolist())
        high_res = filtered.loc[is_high_res]
        # 3. Get the n alleles with the highest dosages + AB probability
        high_res["DS_AB"] = high_res["DS"] + high_res["AB"]

        return high_res.nlargest(1, columns="DS_AB")

    def sort_and_fill(self):
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
                elif len(heterozygous) == 1:
                    # If there is only one allele, look for the second most probable
                    # from the uncertain alleles
                    hetero_low = self._get_ploidy_alleles(
                        alleles, self._ploidy_filter.UNCERTAIN
                    )

                    # from str to output dosage and AB
                    hetero_low_str = hetero_low.apply(
                        lambda row: f"{row.name}({row['DS']}/{row['AB']})", axis=1
                    )

                    to_append = heterozygous.index.tolist() + [hetero_low_str.iloc[0]]
                elif len(heterozygous) == 0:
                    # If there is no heterozygous alleles, look for the two most probable
                    hetero_low = self._get_ploidy_alleles(
                        alleles, self._ploidy_filter.UNCERTAIN, n=2
                    )

                    # from str to output dosage and AB
                    hetero_low_str = hetero_low.apply(
                        lambda row: f"{row.name}({row['DS']}/{row['AB']})", axis=1
                    )

                    to_append = hetero_low_str.tolist()
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
