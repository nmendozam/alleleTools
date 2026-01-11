import argparse
import os
import pytest

from ..argtypes import output_path


class TestOutputPath:
    def test_file_already_exists(self, tmp_path):
        # create a file and pass its path -> should raise because file exists
        existing = tmp_path / "already.alt"
        existing.write_text("x")
        with pytest.raises(argparse.ArgumentTypeError):
            output_path(str(existing))


    def test_ok(self, tmp_path):
        # valid directory, file does not yet exist -> should return the same path
        out = tmp_path / "newdir"
        out.mkdir()
        candidate = out / "out.alt"
        result = output_path(str(candidate))
        assert result == str(candidate)
    
    def test_only_file_in_current_dir(self, tmp_path, monkeypatch):
        # valid current directory, file does not yet exist -> should return the same path
        file = "out.alt"
        monkeypatch.chdir(tmp_path)
        result = output_path(file)
        assert result == str(file)
    
    def test_only_dir_without_file(self, tmp_path):
        # valid directory, file was not provided -> should raise an error
        out = tmp_path / "newdir"
        out.mkdir()
        print(out)
        with pytest.raises(argparse.ArgumentTypeError):
            output_path(str(tmp_path / "newdir") + "/")

    @pytest.mark.parametrize("path", [
            pytest.param(""), # empty string
            pytest.param("."), # current directory only
            pytest.param("/nonexistent_dir/file.alt"), # parent directory does not exist

    ])
    def test_fail(self, path: str):
        with pytest.raises(argparse.ArgumentTypeError):
            output_path(path)
