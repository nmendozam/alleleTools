import pandas as pd


def setup_parser(subparsers):
    parser = subparsers.add_parser(
        name="allele2vcf",
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
        type=str,
        help="path to the file with gene loci information, alternatively you could specify --gene_cluster",
    )
    parser.add_argument(
        "--gene_cluster",
        type=str,
        help="name of the gene cluster (hla or kir), alternatively you could provide a --loci_file"
    )
    parser.add_argument(
        "--vcf",
        type=str,
        help="VCF file to append the alleles to",
        default="file.vcf",
        required=True,
    )
    parser.add_argument(
        "--field_separator",
        type=str,
        help="character separating the fields inside input default is tab",
        default="\t"
    )

    parser.set_defaults(func=call_function)

    return parser

def call_function(args):
    if not (args.gene_cluster or args.loci_file):
        print(
            "Error: either --gene_cluster or --loci_file must be provided."
            "Use -h or --help to see more details."
            )

    genotypes = pd.read_csv(args.input, sep=args.field_separator)

    pairs = _gene_pairs(genotypes.columns)
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
        allele_codes = _diploid_notation(pivotA, pivotB)
        # add column with gene name
        allele_codes["gene"] = pairA
        pre_vcf_alleles = pd.concat([pre_vcf_alleles, allele_codes])

    # Now we add the leading columns and sort the samples (also in columns)
    # to match the base vcf file order.
    # e.g.:
    # CHROM  POS ID  REF  ALT  QUAL  FILTER  INFO  FORMAT  SAMPLE_ID. ...
    gene_loci = pd.read_csv(args.loci_file, sep="\t")
    # column start has chr:pos, divide that into two columns.
    gene_loci[["CHROM", "POS"]] = gene_loci["start"].str.split(":", expand=True)
    # cross gene from pre_vcf_alleles with gene_loci to get the position of each allele.
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

    with open(args.vcf, "a") as f:
        f.write(vcf_alleles.to_csv(index=False, sep="\t", header=False))

def _diploid_notation(pairA, pairB):
    """
    Takes a Multindex series of presence/absence of an allele in a sample
    and returns a single indexed series with the notation for homozygous and
    heterozygous alleles by mixing both tables.
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

def _get_vcf_columns(vcf_file):
    # Read only the line that starts with #CHROM
    with open(vcf_file, "r") as f:
        line = f.readline()
        while not line.startswith("#CHROM"):
            line = f.readline()
    # Remove leading # and \n, then split by tab.
    return line[1:].strip().split("\t")  # [9:]
