# <p align="center"> Gene agnostic AlleleTools! </p>

Working with alleles from highly polymorphic genes, like those in HLA and KIR clusters, is already hard enough! This repo is a collection of tools to facilitate your work on the manipulation and analysis of allele data sets.

The tools in this repo are sorted by category:

1. genotype: a set of tools to facilitate the genotyping process.
2. convert: a group commands that convert allele data between different file formats. Including vcf, csv and our own format .alt (from allele table). You could also convert from vcf to files compatible with [pyHLA](https://github.com/felixfan/PyHLA) and [PyPop](http://pypop.org/index.html).
3. refactor: some commands to normalize allele resolutions and other useful refactoring.

> [!WARNING]
> This project is currently under development in alpha stage. Some things might break. If that is the case, please do not hesitate to file an issue ([clicking here](https://github.com/nmendozam/alleleTools/issues))

## Getting started

To install this package you can use pip:

```bash
pip install pip@git+https://github.com/nmendozam/alleleTools.git
```

It will install the `altools` command in your current environment. So, to execute a command you need to specify three things `altools [tool_category] [tool_name] [input]`.

```bash
altools convert vcf2allele input.vcf
```

## Usage

### Table of Contents

<details>
    <summary>Click to expand</summary>
    1. [Getting Started](#getting-started)
    2. [Usage](#usage)
        - [Convert Genotype to VCF](#convert-genotype-to-vcf)
            - [Genotype file](#genotype-file)
            - [Gene Location List](#gene-location-list)
            - [Template VCF File](#template-vcf-file)
        - [Convert VCF to Genotype Table](#convert-vcf-to-genotype-table)
        - [Normalizing Allele Resolutions](#normalizing-allele-resolutions)
</details>

### Convert genotype to vcf

To convert a genotype file to vcf, you can use the command `allele2vcf`. It will append the genotyped alleles to a vcf file.

```bash
altools convert allele2vcf resources/hla_diversity.txt --loci_file resources/gene_locations.tsv --vcf file_to_append_to.vcf
```

#### Genotype file

The input format is a tab-separated file, where the first column is the sample name and pairs of columns for each gene. The header gene name convention is "gene" + "gene.1". e.g.

```
"id" "sbgroup" "A" "A.1"
"sample1" "CEPH" "03:01" "02:01"
```

#### Gene location list

Additionally the script requires a list of gene locations. The file should be tab-separated with the following format:

```
gene    start
HFE    6:26087441
HLA-A    6:29942554
```

The first column is the gene name and the second column is (chromosome):(position). This position data can be found in [ensembl](https://www.ensembl.org/index.html) or [UCSC](https://genome.ucsc.edu/). The sample file used in this repo was obtained from a post in [IPD-IMGT](https://www.ebi.ac.uk/ipd/imgt/hla/help/genomics.html)

#### Template VCF file

This is a file with the known SNPs and the header. The script will append the genotyped alleles to this file. The header should contain the gene names in the same format as the gene location list. You need to assure that the file contains ONLY the samples in the genotype file and no more. Otherwise the concatenated alleles won't match the header. To filter the samples you can use `bcftools`:

```bash
cut -d' ' -f1 resources/hla_diversity.txt | tail -n +2 | tr -d '"' |uniq > samples_id.txt
bcftools view --force-samples -S samples_id.txt test.vcf > filtered.vcf
```

### Convert vcf to genotype table

Converting from vcf to a genotype table is also useful. For example when the alleles are imputed the output is a vcf file. To convert it to a genotype table you can use the command `vcf2alleles`. The vcf file should be filtered to contain only the HLA genes. You can use `bcftools` to do this. The script will output a .pyhla file that can be used with [pyHLA](https://github.com/felixfan/PyHLA) and [PyPop](http://pypop.org/index.html). The phenotype file is optional should follow the .phe format of plink files.

```bash
bcftools view --include 'ID~"HLA"' raw_imputed.vcf > only_hla.vcf
altools convert vcf2alleles only_hla.vcf --phe input.phe --out output.pyhla
```

### Normalizing allele resolutions

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
bash src/alleleTools/refactor/allele_resolution.sh one file1.tsv file2.tsv
```

### Consensus allele

This script is used to generate a consensus HLA genotype from the result of many
HLA genotyping algorithms. It will generate a vcf file or an allele table.

```bash
altools genotype consensus --input "IKMB_Reports/*.json" --output "output.txt" --format pyhla
```

The input files follow the format of reports generated by the [ikmb HLA genotyping](https://github.com/ikmb/hla) pipeline. These report files should be in a folder called `IKMB_Reports/` in a json format.
