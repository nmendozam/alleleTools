<div style="padding-left: 10%;" align="center">
<p align="center">
    <img src="./docs/source/_static/logo.svg" alt="AlleleTools Icon" width="260"/>
</p>


  <p align="center">
    <br />
    <a href="https://github.com/nmendozam/alleleTools/issues">Report Bug</a>
    &middot;
    <a href="https://alleletools.readthedocs.io">Explore the docs</a>
  </p>
</div>

# <p align="center"> Gene agnostic AlleleTools! </p>

Working with alleles from highly polymorphic genes, like those in HLA and KIR
clusters, is already hard enough! This repo is a collection of tools to
facilitate your work on the manipulation and analysis of allele data sets.

The commands are sorted in three categories:

1. format: a group commands that convert allele data between different file
   formats. Including vcf, csv, our own format .alt (from allele table) and others.
2. plot: some commands to visualize allele data, specially from external sources.
3. test: perform statistical analysis on the allele dataset.

> [!WARNING]
> This project is currently under development in alpha stage. Some things might
> break. If that is the case, please do not hesitate to file an issue
> ([clicking here](https://github.com/nmendozam/alleleTools/issues))

## Getting started

To install this package you can use pip:

```bash
pip install pip@git+https://github.com/nmendozam/alleleTools.git
```

It will install the `altools` command in your current environment. So, to
execute a command you need to specify three things `altools [tool_category]
[tool_name] [input]`.

```bash
altools format from_vcf input.vcf
```

Other commands can be listed using:

```bash
altools --help
altools format --help
```

You can find more information about the usage in the [documentation](https://alleletools.readthedocs.io).

## Cite
If you find this package useful for your research, please consider citing it. You can cite the current repository. Alternatively, although this is a personal project and has not been published yet, the first mention to the package was made in:

Mendoza-Mejía, N., Kolbe, D., Özer, O., Dose, J., Torres, G. G., Franke, A., Nygaard, M., & Nebel, A. (2025). __HLA-DRB1*15:01 is associated with a reduced likelihood of longevity in northern European men__. Genome Medicine, 17(1), 125. https://doi.org/10.1186/s13073-025-01554-1


