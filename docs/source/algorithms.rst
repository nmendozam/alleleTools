Algorithms
==================
Some complex algorithms implemented in AlleleTools are described below.

.. _consensus_algorithm:
Consensus from multiple genotypes
---------------------------------
Given that several solutions to genotype complex genes are not fully accurate,
it is often useful to combine the results of multiple algorithms to obtain a
more reliable consensus genotype. This is a common approach called `N-modular
redundancy`_, where N independent modules (in this case, genotyping algorithms)
provide outputs that are then compared to reach a consensus. 

.. _N-modular redundancy: https://en.wikipedia.org/wiki/Triple_modular_redundancy

AlleleTools implements a voting system to achieve this. Each algorithm's
genotype is considered a "vote" for the alleles it predicts. The consensus
genotype is determined by selecting the alleles that receive the most votes
across all algorithms. This method helps to mitigate the errors of individual
algorithms and provides a more robust genotype.

Why not just use the most accurate algorithm?
    While some algorithms may be more accurate on average, their performance can
    vary depending on the specific dataset, population, and sequencing
    technology used. For instance, some algorithms always provide a genotype,
    even if the loci was covered poorly (only one read) or the covered area is
    not an exact match to the reference. In such cases, it is imposible to know
    if the genotype is correct with only one algorithm. This is especially true
    in scenarios where the data quality and coverage is low. By combining
    multiple algorithms, we can leverage their strengths and mitigate their
    weaknesses.

Why is this a complex algorithm?
    Although initially this looks like a work for a simple string comparison and
    counting, there are several complexities to consider:

    #. *Allele nomenclature*: Different genotyping tools may use different naming conventions for alleles.

    #. *Partial genotypes*: Some genotyping tools may provide genotypes at a lower resolution than others.

    #. *Ambiguities*: Some genotyping tools may return more than two alleles per gene, indicating uncertainty.

    #. *Missing data*: Some genotyping tools may fail to provide genotypes for certain genes.

How does it work?
    Our algorithm addresses these inconsistencies by employing a `search tree`_
    algorithm with a `trie structure`_ to model the nested structure of HLA allele
    nomenclature. 

    The root of the trie is the gene (e.g., HLA-A), and each subsequent level represents an increasing
    field of resolution. The algorithm operates in three steps:

    1. The process begins with the construction of the tries, where a solution trie is built from the reported alleles of each genotyping tool. Every path from the root to a end node represents a specific allele call (e.g., A*01:01 is represented by a single path A → 01 → 01). 

    .. image:: _static/consensus1.svg
       :align: center
       :scale: 50 %
       :width: 600px
       :height: 400px
       :loading: lazy

    2. These tries are merged into a single unified trie. Each node in this merged trie tracks how many tools support that specific allele call. 

    .. image:: _static/consensus2.svg
       :align: center
       :scale: 50 %
       :width: 600px
       :height: 400px
       :loading: lazy

    3. Finally, the algorithm performs consensus calling, going over the unified trie to identify a pair of high-resolution alleles that meet a user-defined threshold (default 0.6). 

    .. _trie structure: https://en.wikipedia.org/wiki/Trie
    .. _search tree: https://en.wikipedia.org/wiki/Search_tree

