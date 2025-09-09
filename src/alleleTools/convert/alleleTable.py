"""
Allele Table Data Structure Module.

This module defines the AlleleTable class, which provides a standardized
data structure for storing and manipulating allele data along with associated
phenotype and covariate information.
"""

import numpy as np
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

    def load_phenotype(self, phenotype_file: str) -> None:
        """
        Load phenotype information from a file into the AlleleTable.

        Args:
            phenotype_file (str): The path to the phenotype file. It should be
                a whitespace-separated values file with a header.  It should
                contain columns "IID" and "phenotype".
        """
        if len(phenotype_file) == 0:
            self.set_phenotype(
                pd.Series(index=self.alleles.index, data=float(0))
            )
            return

        phe = pd.read_csv(phenotype_file, sep=" ", comment="#")

        # Check that the file was formatted correctly
        if "IID" not in phe.columns or "phenotype" not in phe.columns:
            print(
                f"ERROR: Phenotype file ({phenotype_file}) must contain 'IID' and 'phenotype' columns.")
            print("Check if the file is well formatted and witespace-separated")
            exit(1)

        phe.set_index("IID", inplace=True)

        self.set_phenotype(phe["phenotype"])

    def set_phenotype(self, phenotype: pd.Series) -> None:
        """
        Stores the phenotype series in the AlleleTable

        Args:
            phenotype (pd.Series): A series with the phenotypes for the
                AlleleTable. It should have the samples' IDs as index and the
                phenotype values.
        """
        self.phenotype = phenotype
        self.phenotype.name = "phenotype"
        self.phenotype.index.name = "sample"
        self.__verify_phe_samples__()

    def __verify_phe_samples__(self):
        """
        Check if the phenotype series has the same samples all the samples
        that the allele table has.
        """
        missing_samples = self.alleles.index.difference(self.phenotype.index)

        if len(missing_samples) > 0:
            print("ERROR: Phenotype list does not contain all the samples from alleles.")
            print("%d samples are missing:" % len(missing_samples))
            print(missing_samples.to_list())
            exit(1)

    def to_csv(
            self, filename: str, header: bool = True, population: str = ""
    ):
        """
        Export the allele table to a CSV file.

        Args:
            filename (str): The name of the output CSV file.
            header (bool): Flag to store the file with column names or not
            population (str): Adds an extra column in the position left to
                phenotype with a population name. Currently, only one
                population per allele table is supported.
        """
        df = self.alleles.copy()

        # Convert alleles to string
        df = df.map(str).replace("", np.nan)

        if not self.phenotype.empty:
            # Add phenotype to df, checking that the index matches
            df = df.join(self.phenotype, how="inner")
            # Move phenotype to the first column
            col = df.pop("phenotype")
            df.insert(0, "phenotype", col)

        # add population column at the beginning
        if population:
            df.insert(
                0, "population", population
            )

        # move sample id to column
        df.reset_index(inplace=True)

        df.to_csv(
            filename,
            sep="\t",
            index=False,
            na_rep="NA",
            header=header
        )
