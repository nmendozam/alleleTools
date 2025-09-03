import pandas as pd
import pytest

from .ikmb_report import Gene, Report
from .plot_ikmb_coverage import read_reports_asdf, remove_HLA_prefix


def test_remove_HLA_prefix():
    cov = {"HLA-A": [1, 2], "HLA-B": [3, 4], "C": [5]}
    result = remove_HLA_prefix(cov)
    assert "A" in result
    assert "B" in result
    assert "C" in result
    assert "HLA-A" not in result
    assert result["A"] == [1, 2]


def test_gene_mean_coverage():
    coverage = [{"exon": 1, "mean_cov": 10}, {"exon": 2, "mean_cov": 30}]
    calls = {"HLA-HD": ["A*01:01", "A*02:01"]}
    gene = Gene("A", coverage, calls)
    assert gene.mean_coverage() == 20


class TestGeneConsensus:
    @pytest.fixture
    def coverage(self) -> dict:
        return [{"exon": 1, "mean_cov": 10}, {"exon": 2, "mean_cov": 30}]

    def test_similar_calls(self, coverage):
        calls = {
            "alg1": ["A*01:01", "A*02:01"],
            "alg2": ["A*01:01", "A*02:01"],
        }
        gene = Gene("A", coverage, calls)
        assert gene.get_consensus_call() == ("A*01:01,A*02:01", 2)


def test_gene_asdict():
    coverage = [{"exon": 1, "mean_cov": 10}, {"exon": 2, "mean_cov": 30}]
    calls = {"HLA-HD": ["A*01:01", "A*02:01"]}
    gene = Gene("A", coverage, calls)
    d = gene.asdict()
    assert d["gene"] == "A"
    assert d["coverage"] == 20
    assert type(d["consensus"]) == str
    assert type(d["support"]) == int


def test_report_aslist():
    report_dict = {
        "sample": "S1",
        "calls": {
            "A": {"HLA-HD": ["A*01:01", "A*02:01"]},
            "B": {"HLA-HD": ["B*07:02", "B*08:01"]},
        },
        "coverage": {
            "HLA-A": [{"exon": 1, "mean_cov": 10}, {"exon": 2, "mean_cov": 30}],
            "HLA-B": [{"exon": 1, "mean_cov": 40}, {"exon": 2, "mean_cov": 60}],
        },
    }
    report = Report(report_dict)
    aslist = report.aslist()
    assert isinstance(aslist, list)
    assert aslist[0]["gene"] == "A"
    assert aslist[1]["gene"] == "B"
    assert aslist[0]["sample"] == "S1"
    assert aslist[1]["sample"] == "S1"


def test_read_reports_asdf(tmp_path):
    # Prepare two fake json files
    data1 = {
        "sample": "S1",
        "calls": {"A": {"HLA-HD": ["A*01:01", "A*02:01"]}},
        "coverage": {
            "HLA-A": [{"exon": 1, "mean_cov": 10}, {"exon": 2, "mean_cov": 30}]
        },
    }
    data2 = {
        "sample": "S2",
        "calls": {"B": {"HLA-HD": ["B*07:02", "B*08:01"]}},
        "coverage": {
            "HLA-B": [{"exon": 1, "mean_cov": 40}, {"exon": 2, "mean_cov": 60}]
        },
    }
    file1 = tmp_path / "file1.json"
    file2 = tmp_path / "file2.json"
    file1.write_text(str(data1).replace("'", '"'))
    file2.write_text(str(data2).replace("'", '"'))
    df = read_reports_asdf([str(file1), str(file2)])
    assert isinstance(df, pd.DataFrame)
    assert set(df["sample"]) == {"S1", "S2"}
    assert set(df["gene"]) == {"A", "B"}
    assert "coverage" in df.columns
    assert "coverage" in df.columns
    assert "coverage" in df.columns
