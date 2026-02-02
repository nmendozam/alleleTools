.. AlleleTools documentation master file, created by
   sphinx-quickstart on Wed Jan 28 11:22:24 2026.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

AlleleTools: documentation
=======================================
AlleleTools is a CLI (Command Line Interface) tool for researchers working with allele data from complex genes (e.g., HLA, KIR). It aims to provide a comprehensive suite of functionalities for data conversion, quality control, statistical analysis, and visualization.



CLI modules
--------------------

The commands are sorted in categories, referred to as modules:

1. :doc:`Format<modules/format>`: a group commands that convert allele data between different file
   formats. Including vcf, csv, our own format .alt (from allele table) and others.
2. :doc:`Plot<modules/plot>`: some commands to visualize allele data, specially from external sources.
3. :doc:`Test<modules/test>`: perform statistical analysis on the allele dataset.

API Documentation
--------------------

For contributors and developers, an :doc:`API reference<api/alleleTools>` is available.


.. This section is for the left panel contents table:
    Include files here if you want them to appear in the index.

.. toctree::
   :maxdepth: 1
   :hidden:

   file_types
   algorithms

.. toctree::
   :maxdepth: 1
   :caption: CLI modules
   :hidden:

   modules/format
   modules/plot
   modules/test
    

.. toctree::
   :maxdepth: 1
   :caption: API Documentation
   :hidden:

   api/alleleTools