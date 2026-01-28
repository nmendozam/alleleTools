Plot commands
================
The `plot` subcommands create various quick visualizations to facilitate the analysis. Some of these plots are 
useful after performing association analyses to interpret the results.

pleiotropy
----------------
Generates a pleiotropy report of a given allele across multiple phenotypes.

Input dataset
    You need to download the phewas results from `Phewas HLA <https://phewascatalog.org/phewas/#hla>`__.
    Naturally, the file most include summary statistics of the requested allele.


graph_pathogens
-----------------------
Generates a graph useful to compare the proportions of binding and poorly
binding peptides across different pathogens for a given HLA allele. The data is
source from the IEDB database.

Additionally, a fisher's exact test is performed to determine if the proportions
of binding and poorly binding peptides are significantly different between the
pathogens. If a resulting p-value is below 0.05 and the OR is greater than 1, it
indicates that the provided HLA allele may not recognize peptides from one
pathogen as effectively as other alleles of the same locus. This could suggest a
potential susceptibility to infections caused by that pathogen.

.. note::
    Please consider that not all alleles have sufficient data in the IEDB
    database to generate meaningful plots.


plot_ikmb_coverage
-----------------------
Produces a histogram with the coverage of each HLA gene from the ikmb HLA genotyping pipeline.