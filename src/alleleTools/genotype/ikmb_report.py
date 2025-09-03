import json
from typing import List, Tuple

import pandas as pd

from ..allele import Allele


def remove_HLA_prefix(coverage: dict) -> dict:
    """
    Remove HLA- prefix from coverage
    """
    new_coverage = dict()
    for key, value in coverage.items():
        key = key.replace("HLA-", "")
        new_coverage[key] = value
    return new_coverage


def read_json(path: str) -> dict:
    """
    Reads a json file and returns a dict
    """
    file = open(path)
    d = json.load(file)
    file.close()
    return d


class Gene:
    """
    This object holds the genes parsed from genotyping reports
    """

    def __init__(self, name: str, coverage: List[dict], calls: dict):
        self.name = name
        self.coverage = pd.DataFrame(coverage)
        self.calls = calls

        self.alleles = list()
        for alleles in self.calls.values():
            self.alleles.extend([str(Allele(a, gene=self.name)) for a in alleles])

    def __str__(self) -> str:
        return str(
            {
                "name": self.name,
                "calls": self.calls,
                "mean_cov": self.mean_coverage(),
            }
        )

    def mean_coverage(self) -> float:
        """
        Get the mean coverage of exons from the gene
        """
        mean_cov = self.coverage["mean_cov"].mean()
        return float(mean_cov)

    def asdict(self) -> dict:
        """
        Convert this object into a dictionary
        """
        return {
            "gene": self.name,
            "coverage": self.mean_coverage(),
        }


class Report:
    def __init__(self, report: dict):
        self.sample = report["sample"]

        self.__parse_genes__(calls=report["calls"], coverage=report["coverage"])

    def __parse_genes__(self, calls: dict, coverage: dict):
        """
        Parses the genes from the report and adds them to
        `self.genes`
        """
        coverage = remove_HLA_prefix(coverage)
        self.genes = list()

        for gene in coverage.keys():
            g = self.__parse_gene__(
                name=gene, coverage=coverage[gene], calls=calls[gene]
            )
            self.genes.append(g)

    def __parse_gene__(self, name: str, coverage: List[dict], calls: dict):
        """
        Parses a single gene from the report.
        """
        return Gene(name=name, coverage=coverage, calls=calls)

    def aslist(self) -> List[dict]:
        """
        Convert the report into a list of dictionaries
        """
        result = list()
        for gene in self.genes:
            d = gene.asdict()
            d["sample"] = self.sample
            result.append(d)
        return result
