import sys

import pandas as pd


def diploid_notation(pairA, pairB):
    """
    Takes a Multindex series of presence/absence of an allele in a sample
    and returns a single indexed series with the notation for homozygous and
    heterozygous alleles by mixing both tables.
    """
    # assign a code to each allele combination.
    pairA = pairA.astype(int).rename_axis("id")
    pairB = pairB.astype(int).mul(2).rename_axis("id")
    group = pd.concat([pairA, pairB]).groupby("id").sum()
    return group.replace({0: "0|0", 1: "1|0", 2: "0|1", 3: "1|1"})


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


def get_vcf_columns(vcf_file):
    # Read only the line that starts with #CHROM
    with open(vcf_file, "r") as f:
        line = f.readline()
        while not line.startswith("#CHROM"):
            line = f.readline()
    # Remove leading # and \n, then split by tab.
    return line[1:].strip().split("\t")  # [9:]


if __name__ == "__main__":

    genotype_file = sys.argv[1]  # First argument
    gene_loci_file = sys.argv[2]  # Second argument
    vcf_file = sys.argv[3]  # Third argument

    genotypes = pd.read_csv(genotype_file, sep=" ")

    pairs = gene_pairs(genotypes.columns)
    genes = [item for t in pairs for item in t]

    # By this point we have the vcf file without the leading columns.
    # First we pivot the table to have the samples as columns and the alleles as rows.
    # Then we merge both gene alleles into a single row by applying the diploid_notation function.
    # e.g.:
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
        allele_codes = diploid_notation(pivotA, pivotB)
        # Add gene name to the index
        allele_codes.index = "HLA_" + pair[0] + "_" + allele_codes.index
        allele_codes.index.name = "ID"
        allele_codes = allele_codes.reset_index()
        # add column with gene name
        allele_codes["gene"] = "HLA-" + pairA
        pre_vcf_alleles = pd.concat([pre_vcf_alleles, allele_codes])

    # Now we add the leading columns and sort the samples (also in columns)
    # to match the base vcf file order.
    # e.g.:
    # CHROM  POS ID  REF  ALT  QUAL  FILTER  INFO  FORMAT  SAMPLE_ID. ...
    gene_loci = pd.read_csv(gene_loci_file, sep="\t")
    # column start has chr:pos, divide that into two columns.
    gene_loci[["CHROM", "POS"]] = gene_loci["start"].str.split(":", expand=True)
    # cross gene from pre_vcf_alleles with gene_loci to get the position of each allele.
    vcf_alleles = pre_vcf_alleles.merge(gene_loci, how="left", on="gene")
    vcf_alleles = vcf_alleles.assign(
        REF="A", ALT="T", QUAL=".", FILTER="PASS", INFO=".", FORMAT="GT"
    )
    # Remove gene columns and rename
    vcf_alleles = vcf_alleles.drop(columns=["gene", "start"])

    # Sort columns to match the base vcf file order.
    vcf_col = get_vcf_columns(vcf_file)
    vcf_col = [x for x in vcf_col if x in vcf_alleles.columns]
    vcf_alleles = vcf_alleles[vcf_col]

    with open(vcf_file, "a") as f:
        f.write(vcf_alleles.to_csv(index=False, sep="\t", header=False))
