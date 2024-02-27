import pandas as pd


def diploid_notation(row):
    """
    Takes a Multindex series of presence/absence of an allele in a sample
    and returns a single indexed series with the notation for homozygous and
    heterozygous alleles by mixing both tables.
    """
    row.index = row.index.set_names(["allele", "id"])
    row.name = "value"
    num = row["A"].astype(int) + row["A.1"].astype(int).mul(2)
    return num.replace({0: "0|0", 1: "1|0", 2: "0|1", 3: "1|1"})


if __name__ == "__main__":
    genotypes = pd.read_csv("20140702_hla_diversity.txt", sep=" ")

    # By this point we have the vcf file without the leading columns.
    # First we pivot the table to have the samples as columns and the alleles as rows.
    # Then we merge both gene alleles into a single row by applying the diploid_notation function.
    # e.g.:
    # ID  SAMPLE_ID ...
    pivot = genotypes.pivot(index="A", columns="id", values=["A", "A.1"]).notna()
    pre_vcf_alleles = pivot.apply(diploid_notation, axis=1)
    print(pre_vcf_alleles)

    # Now we add the leading columns and sort the samples (also in columns)
    # to match the base vcf file order.
    # e.g.:
    # CHROM  POS ID  REF  ALT  QUAL  FILTER  INFO  FORMAT  SAMPLE_ID. ...
