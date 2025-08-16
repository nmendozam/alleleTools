import pandas as pd


class VCF:
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
        Reads metadata and data frame of the file
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
        Removes a prefix from allele names in the dataframe.
        Usualy it is the name of the gene HLA or KIR
        """
        self.dataframe["ID"] = self.dataframe["ID"].str.replace(prefix, "")

    def get_format(self):
        """
        Get the format of the info column
        """
        formats = self.dataframe["FORMAT"].str.split(":", expand=True)
        return formats.iloc[0].tolist()

    def samples(self):
        columns = set(self.dataframe.columns)
        sample_columns = columns.difference(self.__static_columns)
        return sample_columns

    def samples_dataframe(self):
        return self.dataframe.loc[:, self.samples()]

    def save(self, path: str):
        self.dataframe = pd.DataFrame()
        with open(path, "w") as file:
            file.write(self.metadata)
        self.dataframe.reset_index(inplace=True)
        self.dataframe.rename(columns={"CHROM": "#CHROM"}, inplace=True)
        self.dataframe.to_csv(file, mode="a", header=True)
