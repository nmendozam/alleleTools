import time

import pandas as pd
import requests
from yaspin import yaspin
from yaspin.spinners import Spinners

# ref: This code is based on the examples found in https://github.com/IEDB/IQ-API-use-cases
base_uri='https://query-api.iedb.org'

organism_iris = {
    "human": "NCBITaxon:9606",
}

disease_iris = {
    "autoimmune": "DOID:417",
    "infectious": "DOID:0050117",
    "cancer": "DOID:162",
    "allergy": "DOID:1205",
}

# funciton to print the CURL command given a request
def print_curl_cmd(req):
    url = req.url
    print("curl -X 'GET' '" + url + "'")

def query_mhc(
        hla_allele,
        disease:str = "",
        d_search:bool = False, # Search intensive for parent disease terms
        min_len:int = 0,
        max_len:int = 0,
        source:str = "",
        host:str = ""
        ) -> pd.DataFrame:

    select_columns = [
       "mhc_allele_name",
       "qualitative_measure",
       "source_organism_name",
       "source_organism_iri",
       "structure_description",
       "disease_names",
       "disease_iris",
       "assay_iris",
    ]

    search_params = {
        # Export Type: Full, all data columns
        "select": ",".join(select_columns),
        "order": "assay_iris",
        # Include Positive and Negative Assays
        # "qualitative_measure": "eq.Negative",
        "qualitative_measure": "in.(Negative,Positive,Positive-Low,Positive-Intermediate,Positive-High)",
        # MHC Restriction Type:
        'mhc_allele_name': 'eq.%s' % hla_allele,
        'offset': 0,
    }

    # Disease Data: 
    if disease in disease_iris:
        if d_search: # look for parent terms
            search_params["disease_iri_search"] = 'cs.{"%s"}' % disease_iris["autoimmune"]
        else:
            search_params["disease_iris"] = 'cs.{"%s}' % disease_iris["autoimmune"]
    
    # Epitope structure
    if min_len > 0 or min_len > 0:
        # Only linear epitopes around 8-12 amino acids
        search_params["structure_type"] = "eq.Linear peptide"
        if min_len > 0 and max_len > 0:
            search_params["and"] = "(linear_sequence_length.gte.%s,linear_sequence_length.lte.%s)" % (min_len, max_len)
        elif min_len > 0:
            search_params["linear_sequence_length.gte"] = min_len
        else:
            search_params["linear_sequence_length.lte"] = max_len

    if source in organism_iris:
        search_params["source_organism_iris"] = 'cs.{%s}' % organism_iris[source]

    if host in organism_iris:
        search_params["host_organism_iris"] = 'cs.{%s}' % organism_iris[host]

    table_name='mhc_search' # Only mhc results
    full_url=base_uri + '/' + table_name

    print("Querying IEDB API, this might take a while!")
    chunk_n = 0
    ret_table = pd.DataFrame()
    with yaspin(Spinners.arc, text="( 0 chunks downloaded)") as sp:
        while True:
            sp.text = "(%02d chunks downloaded) waiting for next chunk" % chunk_n
            # print("Search parameters: %s" % json.dumps(search_params, indent=2))
            result = requests.get(full_url, params=search_params)

            if result.status_code != 200:
                print("Failed on offset(%d) [error code %d]:\n" % (chunk_n, result.status_code), result.json())
                print_curl_cmd(result)
                exit()

            if result.json() == []:
                break

            table_chunk = pd.json_normalize(result.json())
            ret_table = pd.concat([ret_table, table_chunk])

            chunk_n += 1
            sp.text = "(%02d chunks downloaded) 2 sec. delay" % chunk_n
            search_params["offset"] += 10000
            time.sleep(2)

    return ret_table

if __name__ == "__main__":
    df = query_mhc("HLA-DRB1*15:01")
    print(df.head())
