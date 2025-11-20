import os
from pathlib import Path
from urllib.request import urlretrieve


def get_asset_path(file: str) -> str:
    dir = Path(__file__).parent.parent.parent.parent
    return os.path.join(dir, "resources", file)

def download_file(url: str, dest: str) -> None:
    try:
        print(f"Downloading file from {url} to {dest}...")
        urlretrieve(url, dest)
    except Exception as e:
        # cleanup partial file if any
        try:
            if os.path.exists(dest):
                os.remove(dest)
        except Exception:
            pass
        raise FileNotFoundError(f"Could not download file from {url}: {e}")