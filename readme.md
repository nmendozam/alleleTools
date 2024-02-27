This script is used to parse typing data and convert append it to a vcf file with snp data. This might be useful for calculating linkage disequilibrium between HLA alleles and SNPs in tagSNP selection and construction of imputation panels.

# Getting started

You can use conda to install the environment to be on se safe side.

```bash
conda env create -f environment.yml
```

But in theory the only requirements are:

- Python 3.6 or higher
- pandas

The sample file with HLA typing is taken from DOI: 10.1371/journal.pone.0097282. To run this script you need to download the file by executing the following command:

```bash
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20140725_hla_genotypes/20140702_hla_diversity.txt
```

> [!NOTE]
> The file requires reformating in the column headers as explained in the input format section.

# Input format

## Typing file

The input format is a tab-separated file, where the first column is the sample name and pairs of columns for each gene. The header gene name convention is "gene" + "gene.1". e.g.

```
"id" "sbgroup" "HLA-A" "HLA-A.1"
"sample1" "CEPH" "03:01" "02:01"
```

## Gene location list

Additionally the script requires a list of gene locations. The file should be a tab-separated file formatted as follows:

```
Gene	Start
HFE	6:26087441
HLA-A	6:29942554
```

Where the first column is the gene name and the second column is (chromosome):(position). This position data can be found in [ensembl](https://www.ensembl.org/index.html) or [UCSC](https://genome.ucsc.edu/). However the sample file used in this repo was obtained from a post in [IPD-IMGT](https://www.ebi.ac.uk/ipd/imgt/hla/help/genomics.html)
