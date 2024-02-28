This script is used to parse genotyping data and append it to a vcf file with snp data. This might be useful for calculating linkage disequilibrium between alleles and SNPs, tagSNP selection, and construction of imputation panels.

# Getting started

You can use conda to install the environment to be on se safe side.

```bash
git clone https://github.com/nmendozam/alleles2vcf.git && cd alleles2vcf
conda env create -f environment.yml
```

But in theory the only requirements are:

- Python 3.6 or higher
- pandas

# Usage

The script requires three input files:

```bash
conda activate vcf
python main.py resources/hla_diversity.txt resources/gene_table.tsv resources/filtered.vcf
```

# Input format

## Typing file

The input format is a tab-separated file, where the first column is the sample name and pairs of columns for each gene. The header gene name convention is "gene" + "gene.1". e.g.

```
"id" "sbgroup" "A" "A.1"
"sample1" "CEPH" "03:01" "02:01"
```

## Gene location list

Additionally the script requires a list of gene locations. The file should be tab-separated with the following format:

```
gene    start
HFE    6:26087441
HLA-A    6:29942554
```

The first column is the gene name and the second column is (chromosome):(position). This position data can be found in [ensembl](https://www.ensembl.org/index.html) or [UCSC](https://genome.ucsc.edu/). The sample file used in this repo was obtained from a post in [IPD-IMGT](https://www.ebi.ac.uk/ipd/imgt/hla/help/genomics.html)

## VCF file

It is required to filter the VCF file, so that it only contains the samples that were genotyped. Otherwise the concatenated alleles won't match the header.

```bash
cut -d' ' -f1 20140702_hla_diversity.txt | tail -n +2 | tr -d '"' |uniq > samples_id.txt
bcftools view --force-samples -S samples_id.txt test.vcf > filtered.vcf
```
