import os 
from .. import alleleTools

def get_asset_path(file:str) -> str:
    path = os.path.dirname(os.path.abspath(alleleTools.__file__))
    # remove two levels to get to the root of the package
    path = os.path.dirname(os.path.dirname(path))

    return os.path.join(path, "resources", file)

