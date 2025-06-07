import argparse
import os

def path(in_string:str):
    if not os.path.exists(in_string):
        raise argparse.ArgumentTypeError(f"{path} does not exist.")

    return in_string


def file_path(in_string:str):
    file = path(in_string) 
    if not os.path.isfile(file):
        raise argparse.ArgumentTypeError(f"{path} is not a valid file")

    return file

