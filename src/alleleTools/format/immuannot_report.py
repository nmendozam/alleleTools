import glob
import gzip
import multiprocessing
import os
import threading
from queue import Queue
from typing import List

import pandas as pd
import progressbar

from ..argtypes import file_path, output_path
from .alleleTable import AlleleTable


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
        help="Convert immuannot reports to allele table format",
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
        help="GTF genotyping files from immuannot. It can be a glob pattern.",
    )
    parser.add_argument(
        "--threads",
        type=int,
        help="Number of parallel process"
    )
    parser.add_argument(
        "--max_template_dist",
        type=int,
        help="Threshold to include alleles on the output table. Calls with high distances are less accurate.",
        default=20,
    )
    parser.add_argument(
        "--phe",
        type=file_path,
        help="input phe file name (to add phenotype column)",
        default="",
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

    def get_attributes(self, feature: str="gene") -> list:
        filtered_df = self.dataframe[self.dataframe.feature == feature]
        raw_attris = filtered_df["attribute"].to_list()

        attris = list()
        for raw in raw_attris:
            attris.append(self.parse_attribute(raw))

        return attris


def get_file_list(input: List[str]) -> List[str]:
    file_list = input

    # Check that input is not a glob
    if len(input) == 1 and not os.path.exists(input[0]):
        file_list = glob.glob(input[0])

    print("Processing %d files" % len(file_list))
    return file_list


def process_gtf_file(file: str, out_queue: Queue) -> None:
    gtf = GTF(file)
    attributes = gtf.get_attributes()

    for item in attributes:
        item["name"] = gtf.name

    out_queue.put(attributes)


def worker(files: Queue, out_queue: Queue) -> None:
    while True:
        file = files.get()
        if file is None:
            files.task_done()
            break

        # Process input file
        process_gtf_file(file, out_queue)

        files.task_done()


def queue_files(files: List[str]) -> Queue:
    tasks = Queue()
    for file in files:
        tasks.put(file)
    return tasks


def print_progress_bar(queue: Queue, max: int) -> None:

    last: float = 0
    b = progressbar.ProgressBar(max_value=max)
    b.start()
    while not queue.empty():
        done = max - queue.qsize()
        if done - last >= max / 1000:
            b.update(done)
            last = done
    b.finish()


def get_results_as_df(out_queue: Queue) -> pd.DataFrame:
    results = list()
    while not out_queue.empty():
        result = out_queue.get()
        results.extend(result)

    return pd.DataFrame(results)

def read_input_files(input_files: List[str], num_workers: int) -> Queue:
    file_list = get_file_list(input_files)
    in_queue = queue_files(file_list)

    out_queue = Queue()

    # Deploy input workers
    threads = []
    for _ in range(num_workers):
        in_queue.put(None)
        thread = threading.Thread(target=worker, args=(in_queue, out_queue))
        thread.start()
        threads.append(thread)

    print_progress_bar(in_queue, max=len(file_list))

    # Wait for all threads to finish
    in_queue.join()
    for thread in threads:
        thread.join()
    
    return out_queue


def call_function(args):
    # Get number of threads
    num_workers = multiprocessing.cpu_count()
    if hasattr(args, "threads") and args.threads:
        num_workers = args.threads
    print("Deploying %d workers" % num_workers)

    # Read input files
    out_queue = read_input_files(args.input, num_workers)

    # Concatenate all samples into a single file
    all = get_results_as_df(out_queue)

    # Filter by max_template_dist if provided
    if args.max_template_dist is not None:
        all = all[all["template_distance"].astype(int) <= args.max_template_dist]

    all[["sample", "strand"]] = all["name"].str.split('.', expand=True)
    all["gene_name"] = all["gene_name"] + '_' + all["strand"]
    table = pd.pivot_table(all, values="template_allele",
                           index="sample", columns="gene_name", aggfunc='sum')
    
    table.set_index(table.index.astype(str), inplace=True)

    alt = AlleleTable()
    alt.alleles = table
    alt.load_phenotype(args.phe)

    alt.to_csv(args.output)
