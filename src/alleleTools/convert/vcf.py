"""
VCF File Handling Module.

This module provides a VCF (Variant Call Format) class for reading, parsing,
and manipulating VCF files containing genetic variant data, particularly
optimized for HLA and KIR allele data.
"""

import pandas as pd


class VCF:
    """
    A class for handling VCF (Variant Call Format) files.

    This class provides methods to read, parse, and manipulate VCF files,
    with specific functionality for handling allele data from polymorphic
    genes like HLA and KIR.

    Attributes:
        metadata (str): VCF header metadata lines
        dataframe (pd.DataFrame): Main VCF data with ID as index

    Args:
        path (str): Path to the VCF file to read
    """

    def __init__(self, path):
        self.metadata = str()
        self.__read_file(path)
        self.__static_columns = [
            "CHROM",
            "POS",
            "REF",
            "ALT",
            "QUAL",
            "FILTER",
            "INFO",
            "FORMAT",
        ]

    def __read_file(self, path):
        """
        Read and parse a VCF file into metadata and dataframe components.

        Separates the VCF header (metadata) from the data section and loads
        the data into a pandas DataFrame with appropriate data types.

        Args:
            path (str): Path to the VCF file to read

        Raises:
            FileNotFoundError: If the specified file doesn't exist
            pd.errors.ParserError: If the file format is invalid
        """
        last_pos = 0
        with open(path, "r") as file:
            # Store header in metadata
            while True:
                line = file.readline()
                if line.startswith("#CHROM"):
                    break
                self.metadata += line
                last_pos = file.tell()

            # Read the rest of the file as df
            file.seek(last_pos)
            self.dataframe = pd.read_csv(
                file,
                sep="\t",
                on_bad_lines="warn",
                dtype={
                    "#CHROM": str,
                    "POS": int,
                    "ID": str,
                    "REF": str,
                    "ALT": str,
                    "QUAL": str,
                    "FILTER": str,
                    "INFO": str,
                    "FORMAT": str,
                },
            )
            self.dataframe.rename(columns={"#CHROM": "CHROM"}, inplace=True)
            self.dataframe.set_index("ID", inplace=True)

    def remove_id_prefix(self, prefix: str):
        """
        Remove a prefix from allele IDs in the dataframe.

        This is commonly used to remove gene prefixes like "HLA_" or "KIR"
        from allele identifiers to standardize naming.

        Args:
            prefix (str): The prefix string to remove from allele IDs

        Example:
            >>> vcf.remove_id_prefix("HLA_")
            # "HLA_A*01:01" becomes "A*01:01"
        """
        self.dataframe.index = self.dataframe.index.str.replace(
            prefix, "", regex=False)

    def get_format(self):
        """
        Extract the format field specification from the VCF.

        Parses the FORMAT column to determine the structure of genotype
        information fields (e.g., GT:DS:AA:AB:BB).

        Returns:
            List[str]: List of format field names in order

        Example:
            >>> vcf.get_format()
            ['GT', 'DS', 'AA', 'AB', 'BB']
        """
        formats = self.dataframe["FORMAT"].str.split(":", expand=True)
        return formats.iloc[0].tolist()

    def samples(self):
        """
        Get the list of sample column names from the VCF.

        Returns all column names that are not part of the standard VCF
        format (i.e., sample-specific genotype columns).

        Returns:
            set: Set of sample column names
        """
        columns = set(self.dataframe.columns)
        sample_columns = columns.difference(self.__static_columns)
        return sample_columns

    def samples_dataframe(self):
        """
        Get a dataframe containing only the sample genotype columns.

        Returns:
            pd.DataFrame: DataFrame with only sample columns, indexed by variant ID
        """
        return self.dataframe.loc[:, self.samples()]

    def save(self, path: str):
        """
        Save the VCF data to a file.

        Writes the metadata header followed by the dataframe in standard
        VCF format.

        Args:
            path (str): Output file path

        Note:
            This method modifies the internal dataframe structure during saving.
        """
        with open(path, "w") as file:
            file.write(self.metadata)

        # Prepare dataframe for output
        output_df = self.dataframe.reset_index()
        output_df.rename(columns={"CHROM": "#CHROM"}, inplace=True)
        output_df.to_csv(path, mode="a", sep="\t", index=False)
