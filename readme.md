This repo has a collection of scripts to convert files from genotyping to vcf and back.

- `allele2vcf.py` parses .tsv files with genotyping data and appends it to a vcf file.
- `vcf2alleles.py` does the opposite, it extracts the genotyping data from a vcf file and writes it to a .tsv file.
- `allele_resolution.sh` normalizes alleles to a uniform resolution of the input file (.tsv).

The output of the latter is compatible with [pyHLA](https://github.com/felixfan/PyHLA) and [PyPop](http://pypop.org/index.html). While the output of the former might be useful for calculating linkage disequilibrium between alleles and SNPs, tagSNP selection, and construction of imputation panels compatible with VCF files.

<!-- TABLE OF CONTENTS -->
<details>
    <summary style="font-size: 24px; font-weight: bold;">Table of Contents</summary>
    <!-- Your content here -->
    <ol>
    <li><a href="#getting-started">Getting Started</a></li>
    <li>
      <a href="#usage">Usage</a>
      <ul>
        <li><a href="#convert-genotype-to-vcf">Convert Genotype to VCF</a></li>
        <ul>
            <li><a href="#genotype-file">Genotype file</a></li>
            <li><a href="#gene-location-list">Gene Location List</a></li>
            <li><a href="#template-vcf-file">Template VCF File</a></li>
        </ul>
        <li><a href="#convert-vcf-to-genotype-table">Convert VCF to Genotype Table</a></li>
        <li><a href="#normalizing-allele-resolutions">Normalizing Allele Resolutions</a></li>
      </ul>
    </li>
  </ol>
</details>

# Getting started

You can use conda to install the environment to be on se safe side.

```bash
git clone https://github.com/nmendozam/alleleTools.git && cd alleleTools
conda env create -f environment.yml
```

But in theory the only requirements are:

- Python 3.6 or higher
- pandas 2.0.3

# Usage

## Convert genotype to vcf

To convert a genotype file to vcf, you can use the script `allele2vcf.py`. It will append the genotyped alleles to the vcf file.

```bash
conda activate vcf
python allele2vcf.py resources/hla_diversity.txt resources/gene_table.tsv resources/template.vcf
```

### Genotype file

The input format is a tab-separated file, where the first column is the sample name and pairs of columns for each gene. The header gene name convention is "gene" + "gene.1". e.g.

```
"id" "sbgroup" "A" "A.1"
"sample1" "CEPH" "03:01" "02:01"
```

### Gene location list

Additionally the script requires a list of gene locations. The file should be tab-separated with the following format:

```
gene    start
HFE    6:26087441
HLA-A    6:29942554
```

The first column is the gene name and the second column is (chromosome):(position). This position data can be found in [ensembl](https://www.ensembl.org/index.html) or [UCSC](https://genome.ucsc.edu/). The sample file used in this repo was obtained from a post in [IPD-IMGT](https://www.ebi.ac.uk/ipd/imgt/hla/help/genomics.html)

### Template VCF file

This is a file with the known SNPs and the header. The script will append the genotyped alleles to this file. The header should contain the gene names in the same format as the gene location list. You need to assure that the file contains ONLY the samples in the genotype file and no more. Otherwise the concatenated alleles won't match the header. To filter the samples you can use `bcftools`:

```bash
cut -d' ' -f1 resources/hla_diversity.txt | tail -n +2 | tr -d '"' |uniq > samples_id.txt
bcftools view --force-samples -S samples_id.txt test.vcf > filtered.vcf
```

## Convert vcf to genotype table

Converting from vcf to a genotype table is also useful. For example when the alleles are imputed the output is a vcf file. To convert it to a genotype table you can use the script `vcf2alleles.py`. The vcf file should be filtered to contain only the HLA genes. You can use `bcftools` to do this. The script will output a .pyhla file that can be used with [pyHLA](https://github.com/felixfan/PyHLA) and [PyPop](http://pypop.org/index.html). The phenotype file is optional should follow the .phe format of plink files.

```bash
conda activate vcf
bcftools view --include 'ID~"HLA"' raw_imputed.vcf > only_hla.vcf
python vcf2alleles.py only_hla.vcf --phe input.phe --out output.pyhla
```

## Normalizing allele resolutions

This script normalizes allele resolutions to a uniform level of the input file, facilitating association analyses. It ensures that alleles, such as 01 and 01:01, which are essentially identical, are recognized as equal by renaming them for consistency.

```
- resolution 1:
    - 01:01 -> 01
    - 01 -> 01
    - 02:03 -> 02
- resolution 2:
    - 01:01 -> 01:01
    - 01 -> 01:01
    - 02:03 -> 02:03
- resolution 3:
    - 01:01 -> 01:01:01
    - 01 -> 01:01:01
    - 02:03 -> 02:03:01
```

If an output file name is not provided, it will be named `*.[resolution]fields.tsv`, containing the resolved alleles. Up to three resolution levels are supported (one, two and three).

```bash
./allele_resolution.sh one file1.tsv file2.tsv
```
