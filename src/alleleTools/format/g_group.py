import math
from typing import Union

import pandas as pd

from ..argtypes import add_out_altable_args
from ..utils.assets import get_asset_path
from .alleleTable import AlleleTable


def setup_parser(subparsers):
    parser = subparsers.add_parser(
        name="g_group",
        description="Group alleles according to the g-group nomenclature",
        epilog="Author: Nicolás Mendoza Mejía (2025)",
    )
    # Input/output arguments
    parser.add_argument(
        "input",
        type=str,
        help="Allele table with the alleles that you want to convert",
    )
    parser = add_out_altable_args(parser)

    parser.set_defaults(func=call_function)

    return parser


def call_function(args):
    group = GrouperHLA()

    df = pd.read_csv(args.input, sep='\t', index_col=0, dtype=str, header=None)

    # pop the first column as phenotype
    phenotype = df.pop(df.columns[0])

    # Melt and apply g-group conversion
    df = df.melt(ignore_index=False, var_name='col_index', value_name='allele')
    df['gene'] = df['allele'].apply(
        lambda x: x.split('*')[0] if pd.notna(x) else x)
    df['g_group'] = df.apply(lambda row: group.lookup(row['gene'], row['allele']), axis=1)

    df.loc[df['g_group'].isna(), 'g_group'] = df.loc[df['g_group'].isna(), 'allele'] 


    alt = AlleleTable()
    alt.alleles = df.pivot_table(index=df.index, columns='col_index', values='g_group', aggfunc='first')
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
            reference_file = get_asset_path("hla_nom_g.txt")

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


    def lookup(self, gene: str, allele: str) -> Union[str, float]:
        if not gene or not allele:
            return float('nan')
        if pd.isna(gene) or pd.isna(allele):
            return float('nan')

        gene_ref = self.index[gene]

        allele = allele.replace(gene + '*', '') # remove gene prefix

        if allele in gene_ref: # return exact match
            ret = gene_ref[allele]
            return gene + '*' + ret if isinstance(ret, str) else ret

        # If not exact match, look for partial matches
        partial = [k for k in gene_ref.keys() if k.startswith(allele)]

        if len(partial) == 0:
            return float('nan')
        
        # If more than one partial match, check if all fit in the same group
        groups = [gene_ref[p] for p in partial]
        uniq_groups = set(groups)

        if len(uniq_groups) == 1: # return single partial match
            ret = uniq_groups.pop()
            return gene + '*' + ret if isinstance(ret, str) else ret
        
        if len(allele.split(':')) > 2: 
            print(f"WARNING: Ambiguous g-group assignment for [{gene}*{allele}]. Multiple possible groups found: {groups}.")

        return float('nan')

