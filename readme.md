This script is used to parse typing data and convert append it to a vcf file with snp data. This might be useful for calculating linkage disequilibrium between HLA alleles and SNPs in tagSNP selection and construction of imputation panels.

# Getting started

The requirements for this script are:

- Python 3.6 or higher
- pandas

The sample file with HLA typing is taken from DOI: 10.1371/journal.pone.0097282. To run this script you need to download the file by executing the following command:

```bash
wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20140725_hla_genotypes/20140702_hla_diversity.txt
```
