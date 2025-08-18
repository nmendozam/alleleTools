import matplotlib.pyplot as plt
import pandas as pd

from ..argtypes import path

from .iedb.epitope_query import organism_iris, query_mhc
from .iedb.taxon_query import query_taxon_ids


def setup_parser(subparsers):
    parser = subparsers.add_parser(
        name='graph_pathogens',
        help='Graph epitope assay results for HLA alleles',
        description='Consult IEDB for assay results on the desired HLA allele. This command produces graphs with the assay results.',
        epilog="Author: Nicolás Mendoza Mejía (2025)",
    )
    parser.add_argument(
        'email',
        type=str,
        help='email required to get the tax ids from NCBI'
    )
    parser.add_argument(
        'allele',
        type=str,
        help='HLA allele that you want to get the data from'
    )
    parser.add_argument(
        "--min_epitope_len",
        type=int,
        help="Minimum epitope length to query from the database",
        default=0
    )
    parser.add_argument(
        "--max_epitope_len",
        type=int,
        help="Maximum epitope length to query from the database",
        default=0
    )
    parser.add_argument(
        "--source",
        type=int,
        help=f"Epitopes' source organism. Supported organisms are: [{organism_iris.keys()}]. If your organism of interest is not supported, please submit a request in issues. It is a quick change in the code.",
        default=0
    )
    parser.add_argument(
        "--host",
        type=int,
        help=f"Epitopes' host organism. Supported organisms are: [{organism_iris.keys()}]. If your organism of interest is not supported, please submit a request in issues. It is a quick change in the code.",
        default=0
    )
    parser.add_argument(
        "--output_basename",
        type=path,
        help="Base name of the output files, since this command outputs multiple graphs.",
        default="output_graphs"
    )

    parser.set_defaults(func=call_function)

    return parser


def call_function(args):
    # Query the IEDB API for a specific HLA allele
    data = query_mhc(
        args.allele,
        min_len=args.min_epitope_len,
        max_len=args.max_epitope_len
    )

    # Check which rows contain the NCBI IDs
    contains_ncbi_id = data["source_organism_iri"].str.contains(
        "NCBITaxon", case=False, na=False)

    # Pre-process taxa id
    leading_id = "NCBITaxon:"
    data["TaxId"] = data["source_organism_iri"].replace(
        leading_id, "", regex=True)
    taxon_ids = data.loc[contains_ncbi_id,
                         "TaxId"].dropna().unique().astype(int)
    taxon_ids = map(str, taxon_ids)

    # Proceed to get the genus and family of each ncbi taxon id
    taxon_ranks = query_taxon_ids(taxon_ids, args.email)

    merged = data.merge(taxon_ranks, left_on="TaxId", right_on="TaxId")

    graph_by_genus(merged, "Bacteria",
                   f"{args.output_basename}_bacteria_genus.svg")
    graph_by_genus(merged, "Viruses",
                   f"{args.output_basename}_virus_genus.svg")

# %% Stacked plot by percentage


def adjustFigAspect(fig, aspect=1):
    """
    Adjust the subplot parameters so that the figure has the correct
    aspect ratio.
    """
    xsize, ysize = fig.get_size_inches()
    minsize = min(xsize, ysize)
    xlim = 0.4 * minsize / xsize
    ylim = 0.4 * minsize / ysize
    if aspect < 1:
        xlim *= aspect
    else:
        ylim /= aspect
    fig.subplots_adjust(
        left=0.5 - xlim, right=0.5 + xlim, bottom=0.5 - ylim, top=0.5 + ylim
    )


def graph_by_genus(data: pd.DataFrame, division: str, output_file: str):

    grouped = (
        data[data["Division"] == division]
        .groupby(["genus", "qualitative_measure"])
        .count()
    )
    grouped = grouped["assay_iris"].unstack()
    n_samples = grouped.sum(axis=1)
    grouped.index = [
        "{} (n={})".format(y.replace("_", " ").title(), int(val))
        for y, val in n_samples.to_dict().items()
    ]
    grouped = grouped[grouped.sum(axis=1) > 10]

    # sort by the sum of the columns
    grouped = grouped.loc[grouped.sum(
        axis=1).sort_values(ascending=False).index]

    # normalize rows
    grouped = grouped.div(grouped.sum(axis=1), axis=0)

    # Filter qualitative values
    # expected_values = set([
    #         "Positive-Low",
    #         "Negative",
    #         "Positive",
    #         "Positive-Intermediate",
    #         "Positive-High",
    #     ]
    # )
    # qual_values = grouped.index.get_level_values("qualitative_measure")
    # qual_values = expected_values.intersection(set(qual_values))
    # grouped = grouped.loc[:, qual_values]

    # grouped = grouped.sort_index(axis=1)

    grouped[["Negative", "Positive-Low"]
            ] = grouped[["Negative", "Positive-Low"]] * -1

    colors = ["#8b0000", "#f6655f", "#fdb966", "#97e692", "#069d59"]
    fig, ax = plt.subplots()
    adjustFigAspect(fig, aspect=0.7)
    grouped.plot(kind="barh", stacked=True, width=0.7, color=colors, ax=ax)
    plt.legend(loc="upper left", bbox_to_anchor=(1, 1))
    # set the x-axis limits to -1, 1
    ax.set_xlim(-1, 1)
    ax.autoscale(enable=True, axis="y", tight=True)

    plt.savefig(output_file, bbox_inches="tight")
