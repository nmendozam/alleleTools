"""
Given that we want to group epitopes by genus, we need to
query the Taxonomy database to get the genus name for each
TaxId.
"""

from typing import List, cast

import pandas as pd
from Bio import Entrez
from Bio.Entrez import efetch, read


def extract_rank(row) -> pd.Series:
    linage = row["LineageEx"]
    extract = ["TaxId", "ScientificName", "Division"]
    new_row = {tag: row[tag] for tag in extract}
    for taxon in linage:
        rank = taxon["Rank"]
        new_row[rank] = taxon["ScientificName"]
    return pd.Series(new_row)


def query_taxon_ids(taxon_ids: List[str], email: str) -> pd.DataFrame:
    """
    Given a list of taxon ids, query the NCBI taxonomy database
    and return a DataFrame with the results.
    """
    Entrez.email = email

    taxon_id = ",".join(taxon_ids)
    handle = efetch(db="taxonomy", id=str(taxon_id), retmode="xml")
    records = read(handle)
    df = pd.json_normalize(cast(List[dict], records))

    return df.apply(extract_rank, axis=1)


if __name__ == "__main__":
    taxon_ids = ["1773"]

    email = ""

    df = query_taxon_ids(taxon_ids, email)
    print(df)
