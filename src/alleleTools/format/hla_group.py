import os
from typing import Union

import pandas as pd

from ..argtypes import add_out_altable_args
from ..utils.assets import download_file, get_asset_path
from ..allele import AlleleParser
from .alleleTable import AlleleTable
import swifter


def setup_parser(subparsers):
    parser = subparsers.add_parser(
        name="hla_group",
        help="Group alleles by g- or p-group nomenclature",
        description="Group alleles according to the g- or p-group nomenclature",
        epilog="Author: Nicolás Mendoza Mejía (2025)",
    )
    # Input/output arguments
    parser.add_argument(
        "input",
        type=str,
        help="Allele table with the alleles that you want to convert",
    )
    parser.add_argument(
        "--group_type",
        type=str,
        choices=["g-group", "p-group"],
        help="Type of group nomenclature to use (g-group or p-group)",
        default="g-group",
    )
    parser = add_out_altable_args(parser)

    parser.set_defaults(func=call_function)

    return parser


def call_function(args):
    group = GrouperHLA(args.group_type)

    alt_before = AlleleTable.open(args.input)

    # Melt and apply group conversion
    df = alt_before.alleles.melt(ignore_index=False, var_name='col_index', value_name='allele')
    df.index.name = 'sample'
    df = df.reset_index()

    # Apply conversions
    df[['gene', 'allele_stripped']] = df['allele'].str.split('*', n=1, expand=True)
    df = quick_lookup(group, df)

    # Slow lookup on remaining NaNs for alleles with a resolution 2 or more fields
    miss_group = df['group'].isna() & (df['allele_stripped'].str.count(':') >= 1)
    missing = df[miss_group] 
    df.loc[miss_group, 'group'] = missing.swifter.apply(lambda row: group.lookup(row['gene'], row['allele_stripped']), axis=1)

    df.loc[df['group'].isna(), 'group'] = df.loc[df['group'].isna(), 'allele'] 


    alt = AlleleTable()
    alt.alleles = df.pivot_table(index="sample", columns='col_index', values='group', aggfunc='first')
    alt.phenotype = alt_before.phenotype
    alt.phenotype.name = "phenotype"
    alt.to_csv(args.output)

def quick_lookup(group: "GrouperHLA", df: pd.DataFrame) -> pd.DataFrame:
    """
    This optimization was not comprehensibly tested to see if it improves performance, but
    in theory it should be faster than doing individual lookups. It only does exact matches though.

    Args:
        group (GrouperHLA): The GrouperHLA instance with the g-group/p-group index.
        df (pd.DataFrame): DataFrame with 'gene' and 'allele_stripped' columns.
    Returns:
        pd.DataFrame: DataFrame with an added 'group' column.
    """
    ref_rows = []
    for g, mapping in group.index.items():          # mapping: {allele_stripped: group_value}
        for allele_s, grp in mapping.items():
            ref_rows.append((g, allele_s, grp))
    ref_df = pd.DataFrame(ref_rows, columns=['gene', 'allele_stripped', 'group'])

    # Bulk join exact matches (fast)
    return df.merge(ref_df, how='left', on=['gene', 'allele_stripped'])

class GrouperHLA:
    def __init__(self, reference_file: str = "g-group"):
        """
        Load the reference file that was downloaded from:
        https://hla.alleles.org/pages/wmda/g_groups/
        """
        file_path = self._get_group_norm_file(reference_file)

        ref = pd.read_csv(file_path,
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
                "https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/wmda/hla_nom_g.txt"
                ),
            "p-group":(
                "hla_nom_p.txt",
                "https://raw.githubusercontent.com/ANHIG/IMGTHLA/Latest/wmda/hla_nom_p.txt"
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

    def lookup(self, gene: str, allele_stripped: str) -> Union[str, None]:
        """
        Lookup g-group for the given allele. If exact match is not found,
        attempt to find a partial match.
        """
        gene_ref = self.index.get(gene, None)
        if gene_ref is None or allele_stripped is None:
            return None

        if allele_stripped not in gene_ref:
            return self.lookup_partial(gene, allele_stripped)

        return self.lookup_exact(gene, allele_stripped)

    def lookup_exact(self, gene: str, allele_stripped: str) -> Union[str, None]:
        """
        Finds a g-group for the exact allele provided.
        """
        gene_ref = self.index.get(gene, None)
        if gene_ref is None or allele_stripped is None:
            return None

        ret = gene_ref[allele_stripped]
        return gene + '*' + ret if isinstance(ret, str) else ret

    def lookup_partial(self, gene: str, allele_stripped: str) -> Union[str, None]:
        """
        Finds a g-group for the partial allele provided. The index does not
        contain all posible combinations of alleles. It only has the highest
        resolution alleles, so looking for a partial match is necessary.

        A g-group is returned only if one possible match is found.
        """
        # TODO: Optimize this function to avoid looping over all entries
        # could be done by building a trie structure for the alleles.

        gene_ref = self.index.get(gene, None)
        if gene_ref is None or allele_stripped is None:
            return None

        # Look for partial matches
        partial = [k for k in gene_ref.keys() if k.startswith(allele_stripped)]

        if len(partial) == 0:
            return None
        
        # If more than one partial match, check if all fit in the same group
        groups = [gene_ref[p] for p in partial]
        uniq_groups = set(groups)

        if len(uniq_groups) != 1: # Ambiguous match
            return self._resolve_ambiguity(gene, allele_stripped, uniq_groups)

        # Return the single partial match
        ret = uniq_groups.pop()
        return gene + '*' + ret if isinstance(ret, str) else ret
        
    
    def _resolve_ambiguity(self, gene: str, allele_stripped: str, uniq_groups: set) -> Union[str, None]:
        parser = AlleleParser("hla")
        allele_p = parser.parse(gene + '*' + allele_stripped)

        # Alleles of one field or less shouldn't be resolved
        if len(allele_p) <= 1:
            return None

        groups_p = [parser.parse(gene + '*' + group) for group in uniq_groups if isinstance(group, str)]

        truncated_groups = set([p.truncate(len(allele_p)) for p in groups_p])

        if len(truncated_groups) != 1:
            return None
        
        return str(truncated_groups.pop())

