import os
import re
from typing import Tuple, Union

import pandas as pd

from ..argtypes import add_out_altable_args
from ..utils.assets import download_file, get_asset_path
from .alleleTable import AlleleTable


def setup_parser(subparsers):
    parser = subparsers.add_parser(
        name="hla_group",
        description="Group alleles according to the g-group nomenclature",
        epilog="Author: Nicolás Mendoza Mejía (2025)",
    )
    # Input/output arguments
    parser.add_argument(
        "input",
        type=str,
        help="Allele table with the alleles that you want to convert",
    )
    parser.add_argument(
        "group_type",
        type=str,
        choices=["g-group", "p-group"],
        help="Type of group nomenclature to use (g-group or p-group)",
    )
    parser = add_out_altable_args(parser)

    parser.set_defaults(func=call_function)

    return parser


def call_function(args):
    group = GrouperHLA()

    df = pd.read_csv(args.input, sep='\t', index_col=0, dtype=str, header=None)

    # pop the first column as phenotype
    phenotype = df.pop(df.columns[0])

    # Melt and apply group conversion
    df = df.melt(ignore_index=False, var_name='col_index', value_name='allele')
    df['gene'] = df['allele'].apply(
        lambda x: x.split('*')[0] if pd.notna(x) else x)
    df['group'] = df.apply(lambda row: group.lookup(row['gene'], row['allele']), axis=1)

    df.loc[df['group'].isna(), 'group'] = df.loc[df['group'].isna(), 'allele'] 


    alt = AlleleTable()
    alt.alleles = df.pivot_table(index=df.index, columns='col_index', values='group', aggfunc='first')
    alt.phenotype = phenotype
    alt.phenotype.name = "phenotype"
    alt.to_csv(args.output)

class GrouperHLA:
    def __init__(self, reference_file: str = ""):
        """
        Load the reference file that was downloaded from:
        https://hla.alleles.org/pages/wmda/g_groups/
        """
        if not reference_file:
            reference_file = self._get_group_norm_file("g-group")

        ref = pd.read_csv(reference_file,
                        sep=';', comment="#", header=None)
        ref.columns = ["gene", "alleles", "groups"]

        ref["gene"] = ref["gene"].str.replace("*", "")
        ref["alleles"] = ref["alleles"].apply(lambda x: x.split('/'))

        # generate and index for each gene
        index = dict()
        grouped = ref.groupby("gene")

        for gene, group in grouped:
            index[gene] = dict()
            dic = group.apply(
                lambda row: {a: row["groups"] for a in row["alleles"]}, axis=1).to_list()
            dic = {k: v for d in dic for k, v in d.items()}
            index[gene].update(dic)

        self.index = index
    
    def _get_group_norm_file(self, group: str) -> str:
        group_files = {
            "g-group":(
                "hla_nom_g.txt",
                "https://hla.alleles.org/wmda/g_groups/hla_nom_g.txt"
                ),
            "p-group":(
                "hla_nom_p.txt",
                "https://hla.alleles.org/wmda/p_groups/hla_nom_p.txt"
                ),
            }

        if group not in group_files:
            raise ValueError(f"Unknown group type: {group}. Supported types are {list(group_files.keys())}.")

        file, url = group_files[group]
        path = get_asset_path(file)

        if not os.path.exists(path):
            download_file(
                url=url,
                dest=path,
            )
        return path
    
    def _prepare(self, gene: str, allele: str) -> Tuple[Union[dict, None], Union[str, None]]:
        """
        Validate inputs and return (gene_ref, allele_stripped).
        Returns (None, None) if validation fails (caller will return NaN).
        """
        if not gene or not allele:
            return None, None
        if pd.isna(gene) or pd.isna(allele):
            return None, None
        if gene not in self.index:
            return None, None

        gene_ref = self.index[gene]
        allele_stripped = allele.replace(gene + '*', '')

        return gene_ref, allele_stripped

    def lookup(self, gene: str, allele: str) -> Union[str, None]:
        """
        Lookup g-group for the given allele. If exact match is not found,
        attempt to find a partial match.
        """
        gene_ref, allele_stripped = self._prepare(gene, allele)
        if gene_ref is None or allele_stripped is None:
            return None

        if allele_stripped not in gene_ref:
            return self.lookup_partial(gene, allele)

        return self.lookup_exact(gene, allele)

    def lookup_exact(self, gene: str, allele: str) -> Union[str, None]:
        """
        Finds a g-group for the exact allele provided.
        """
        gene_ref, allele_stripped = self._prepare(gene, allele)
        if gene_ref is None or allele_stripped is None:
            return None

        ret = gene_ref[allele_stripped]
        return gene + '*' + ret if isinstance(ret, str) else ret

    def lookup_partial(self, gene: str, allele: str) -> Union[str, None]:
        """
        Finds a g-group for the partial allele provided. The index does not
        contain all posible combinations of alleles. It only has the highest
        resolution alleles, so looking for a partial match is necessary.

        A g-group is returned only if one possible match is found.
        """
        gene_ref, allele_stripped = self._prepare(gene, allele)
        if gene_ref is None or allele_stripped is None:
            return None

        # Look for partial matches
        partial = [k for k in gene_ref.keys() if k.startswith(allele_stripped)]

        if len(partial) == 0:
            return None
        
        # If more than one partial match, check if all fit in the same group
        groups = [gene_ref[p] for p in partial]
        uniq_groups = set(groups)

        if len(uniq_groups) == 1: # return single partial match
            ret = uniq_groups.pop()
            return gene + '*' + ret if isinstance(ret, str) else ret
        
        # 1-field alleles tend to be ambiguous, so just warn for alleles with more fields.
        if len(allele_stripped.split(':')) > 1: 
            print(f"WARNING: Ambiguous g-group assignment for [{gene}*{allele_stripped}]. Multiple possible groups found: {groups}.")

        return None

