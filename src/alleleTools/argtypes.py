import argparse
import os
from pathlib import Path


def path(in_string:str):
    if not os.path.exists(in_string):
        raise argparse.ArgumentTypeError(f"{path} does not exist.")

    return in_string

def output_path(in_string:str):
    basepath = Path(in_string).parent
    if not basepath.exists():
        raise argparse.ArgumentTypeError(f"the directory {basepath} does not exist.")
    elif os.path.isfile(in_string):
        raise argparse.ArgumentTypeError(f"{in_string} already exists!")

    return in_string

def file_path(in_string:str):
    file = path(in_string) 
    if not os.path.isfile(file):
        raise argparse.ArgumentTypeError(f"{path} is not a valid file")

    return file

def csv_file(in_string:str):
    file = file_path(in_string) 
    filename, file_extension = os.path.splitext(file)
    if not file_extension.lower() == '.csv':
        raise argparse.ArgumentTypeError(f"{path} is not a valid CSV file")

    return file

