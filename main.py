import pandas as pd


def diploid_notation(row):
    """
    Takes a Multindex series of presence/absence of an allele in a sample
    and returns a single indexed series with the notation for homozygous and
    heterozygous alleles by mixing both tables.
    """
    row.index = row.index.set_names(["allele", "id"])
    row.name = "value"
    # Get the name of allele pairs from the index.
    pair = row.index.get_level_values("allele").unique().to_list()
    # assign a code to each allele combination.
    num = row[pair[0]].astype(int) + row[pair[1]].astype(int).mul(2)
    return num.replace({0: "0|0", 1: "1|0", 2: "0|1", 3: "1|1"})


def gene_pairs(lst):
    """
    Takes a list and return a list of tuples with the elements in pairs,
    that follow the gene gene.1 convention.
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


if __name__ == "__main__":
    genotypes = pd.read_csv("20140702_hla_diversity.txt", sep=" ")

    pairs = gene_pairs(genotypes.columns)
    genes = [item for t in pairs for item in t]

    # By this point we have the vcf file without the leading columns.
    # First we pivot the table to have the samples as columns and the alleles as rows.
    # Then we merge both gene alleles into a single row by applying the diploid_notation function.
    # e.g.:
    # ID  SAMPLE_ID ...
    pre_vcf_alleles = pd.DataFrame()
    pivot = genotypes.pivot(index="A", columns="id", values=genes).notna()
    for pair in pairs:
        allele_codes = pivot[pair].apply(diploid_notation, axis=1)
        # Add gene name to the index
        allele_codes.index = "HLA_" + pair[0] + "_" + allele_codes.index
    pre_vcf_alleles = pd.concat([pre_vcf_alleles, allele_codes], axis=1)

    # Now we add the leading columns and sort the samples (also in columns)
    # to match the base vcf file order.
    # e.g.:
    # CHROM  POS ID  REF  ALT  QUAL  FILTER  INFO  FORMAT  SAMPLE_ID. ...
