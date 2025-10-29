import os 
import inspect
from pathlib import Path

def get_asset_path(file: str) -> str:
    this_file = inspect.getfile(inspect.currentframe())
    dir = Path(this_file).parent.parent.parent.parent
    return os.path.join(dir, "resources", file)