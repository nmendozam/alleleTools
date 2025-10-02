import glob
import gzip
import os

import pandas as pd

from ..argtypes import output_path


def setup_parser(subparsers):
    """
    Set up the argument parser for reading and converting immuannot
    report.

    Args:
        subparsers: The subparsers object to add this command to.

    Returns:
        argparse.ArgumentParser: The configured parser for consensus.
    """
    parser = subparsers.add_parser(
        name="from_immuannot",
        description="""
        This command converts a group of file reports from immuannot and
        converts the genotyping data to allele table.
        """,
        epilog="Author: Nicolás Mendoza Mejía (2025)",
    )
    parser.add_argument(
        "input",
        metavar="path",
        type=str,
        nargs="+",
        help="GTF genotyping files from immuannot",
    )
    parser.add_argument(
        "--output",
        metavar="path",
        type=output_path,
        help="Path to output file",
        default="output.alt",
    )

    parser.set_defaults(func=call_function)

    return parser


class GTF:
    def __init__(self, path: str):
        self.metadata = str()
        self.dataframe = pd.DataFrame()

        self.name = self.__get_name__(path)
        self.__open_file__(path)

    def __get_name__(self, path: str) -> str:
        file_name = os.path.basename(path)

        name = file_name.replace('.gz', '')
        name = name.replace('.gtf', '')

        return name

    def __open_file__(self, path: str) -> None:
        file = gzip.open(path, 'rt')

        # Store header in metadata
        last_pos = 0
        while True:
            line = file.readline()
            if not line.startswith("##"):
                break
            self.metadata += line
            last_pos = file.tell()
        # Read the rest of the file as data frame
        file.seek(last_pos)
        self.dataframe = pd.read_csv(
            file,
            sep="\t",
            on_bad_lines="warn",
            header=None,
        )

        self.dataframe.columns = [
            "seqname",
            "source",
            "feature",
            "start",
            "end",
            "score",
            "strand",
            "frame",
            "attribute",
        ]
        file.close()

    def parse_attribute(self, line: str) -> dict:
        items = line.split('; ')

        ret = dict()
        for item in items:
            key, value = item.split(' ')
            ret[key] = value.replace('"', '').replace(';', '')

        return ret

    def get_attributes(self, feature="gene") -> list:
        filtered_df = self.dataframe[self.dataframe.feature == feature]
        raw_attris = filtered_df["attribute"].to_list()

        attris = list()
        for raw in raw_attris:
            attris.append(self.parse_attribute(raw))

        return attris


def get_file_list(input):
    if len(input) == 1 and not os.path.exists(input[0]):
        file_list = glob.glob(input[0])
        print("processing from glob: ", len(file_list))
        return file_list
    print("processing: ", len(input))
    return input


def call_function(args):
    all = pd.DataFrame()

    file_list = get_file_list(args.input)
    for file in file_list:
        gtf = GTF(file)
        attributes = gtf.get_attributes()

        df = pd.DataFrame(attributes)

        name, strand = gtf.name.split('.')
        df["sample"] = name
        df["strand"] = strand

        all = pd.concat([all, df])

    all["gene_name"] = all["gene_name"] + '.' + all["strand"]
    table = pd.pivot_table(all, values="template_allele",
                           index="sample", columns="gene_name", aggfunc='sum')

    table.to_csv(args.output, sep='\t')
