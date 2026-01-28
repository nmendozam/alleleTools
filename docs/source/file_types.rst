File formats
========================================
Various file formats are supported. Here we provide an overview of the most common. Specially those that are not described in detail elsewhere in other packages.


.. _alt_format:

allele table (.alt)
------------------------------
The allele table format (.alt) is a simple tab-delimited text file designed to
represent allele data in a structured manner. Each row corresponds to a sample,
and each column represents a specific gene or locus. The first column contains
sample identifiers. The second column contains additional phenotype information.
The subsequent columns contain allele information for each gene. There most be
two columns per gene, representing the two alleles for diploid organisms.
::

    sample    phenotype      HLA-A_1    HLA-A_2
    Sample1     1            A*01:01    A*02:01
    Sample2     0            A*03:01    NA

Missing alleles must be represented as "NA", otherwise they will be treated as
valid alleles depending on the `allele_parsers` used.

.. note::
    Certain commands may parse the alleles to perform specific tasks. Only
    parsers for HLA and KIR are supported by default. Please refer to
    :ref:`allele_parsers` for more information on how to provide custom parsers.

.. _allele_parsers:

allele parsers
------------------------------
AlleleTools includes built-in parsers for HLA and KIR allele nomenclature. These
parsers understand the hierarchical structure of allele names. Additional parsers
can be added by users to support other gene systems.

The config is a json file that specifies the parsing rules for different gene systems.
Two parsing mechanisms are supported regex and delimited. Regardless of the mechanism,
the strategy is to determine the gene name and split the fields that compose the allele.

delimited parser
    This is the simplest parser. However, it only works for gene systems with a
    nomenclature that uses special characters to delimit the gene name and the fields.
    For example, HLA alleles use an asterisk (*) to separate the gene from the fields,
    and colons (:) to separate the fields. In contrast, KIR alleles use astersk to separate
    gene from fields, but do not use any delimiter between fields. In this case,
    the regex must be used. 

    As an example, here is the delimited parser for HLA alleles:
    ::

        {
        "hla_delimited": {
            "type": "delimited",
            "gene_delimiter": "*",
            "field_delimiter": ":"
        }
        }

regex parser
    Given the complexity of regular expressions, we recommend using 
    delimited parsers when possible. However, for gene systems with more complex
    nomenclature, the regex parser provides greater flexibility. Please nota that we
    use the re module from python, so you should follow its `syntax <https://docs.python.org/3/library/re.html#regular-expression-syntax>`__.
    
    This system uses group names to identify the gene and fields. The group names must be:

    - gene: the gene name
    - field1, field2, field3, field4: the fields of the allele
    - suffix: any suffix at the end of the allele name

    As an example, here is the regex parser for HLA alleles:
    ::

        {
        "hla": {
            "type": "regex",
            "pattern":"(?:(?P<gene>\\w+)\\*)(?P<field1>\\d{2})(?::(?P<field2>\\d{2}))?(?::(?P<field3>\\d{2}))?(?::(?P<field4>\\d{2}))?(?P<suffix>\\w)?"
        }
        }