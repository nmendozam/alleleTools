import os 
import inspect
from pathlib import Path

def get_asset_path(file: str) -> str:
    dir = Path(__file__).parent.parent.parent.parent
    return os.path.join(dir, "resources", file)