"""
Allele Table Data Structure Module.

This module defines the AlleleTable class, which provides a standardized
data structure for storing and manipulating allele data along with associated
phenotype and covariate information.
"""

import pandas as pd


class AlleleTable:
    """
    A standardized data structure for storing allele data with metadata.

    This class provides a unified interface for handling allele data from
    polymorphic genes, along with associated phenotype and covariate information.
    It serves as the core data structure for allele analysis workflows.

    Attributes:
        alleles (pd.DataFrame): Main allele data with samples as rows and
           genes/alleles as columns
        phenotype (pd.Series): Phenotype information indexed by sample ID
        covariates (pd.DataFrame): Covariate data with samples as rows and
          covariates as columns

    Example:
        >>> table = AlleleTable()
        >>> # Load allele data
        >>> table.alleles = pd.DataFrame(...)
        >>> # Add phenotype information
        >>> table.phenotype = pd.Series(...)
    """

    def __init__(self):
        """
        Initialize an empty AlleleTable.

        Creates empty pandas structures for alleles, phenotypes, and covariates
        that can be populated with data.
        """
        self.alleles = pd.DataFrame()
        self.phenotype = pd.Series()
        self.covariates = pd.DataFrame()
