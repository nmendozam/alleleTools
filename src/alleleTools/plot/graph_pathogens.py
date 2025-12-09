"""
Pathogen Recognition Visualization Module.

This module provides functionality to query the IEDB (Immune Epitope Database)
for visualizing the results of epitope assays related to a specific HLA allele.
It helps researchers understand which pathogens might the provided allele bind to,
potentially providing resistance against or susceptibility to infections.

Data source: IEDB (Immune Epitope Database) - https://www.iedb.org/

Author: Nicolás Mendoza Mejía (2025)
"""

import matplotlib.pyplot as plt
import pandas as pd

from ..argtypes import path
from .iedb.epitope_query import organism_iris, query_mhc
from .iedb.taxon_query import query_taxon_ids


def setup_parser(subparsers):
    """
    Set up the argument parser for the graph_pathogens command.

    Args:
        subparsers: The subparsers object to add this command to.

    Returns:
        argparse.ArgumentParser: The configured parser for pathogen visualization.
    """
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
    """
    Main function to execute pathogen association analysis and visualization.

    This function orchestrates the complete workflow:
    1. Queries IEDB for epitope-MHC binding data
    2. Extracts and processes NCBI taxonomy IDs
    3. Retrieves taxonomic information from NCBI
    4. Generates separate visualizations for bacteria and viruses

    Args:
        args: Parsed command line arguments containing:
            - email: Email for NCBI API access
            - allele: HLA allele to analyze
            - min_epitope_len: Minimum epitope length filter
            - max_epitope_len: Maximum epitope length filter
            - source: Source organism filter
            - host: Host organism filter
            - output_basename: Base name for output files
    """
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
    data["TaxId"] = data["source_organism_iri"]\
        .replace( leading_id, "", regex=True)
    taxon_ids = data.loc[contains_ncbi_id, "TaxId"]\
        .dropna().unique().astype(int)
    taxon_ids = taxon_ids.astype(str).tolist()

    # Proceed to get the genus and family of each ncbi taxon id
    taxon_ranks = query_taxon_ids(taxon_ids, args.email)

    merged = data.merge(taxon_ranks, left_on="TaxId", right_on="TaxId")

    graph_by_genus(merged, "Bacteria",
                   f"{args.output_basename}_bacteria_genus.svg")
    graph_by_genus(merged, "Viruses",
                   f"{args.output_basename}_virus_genus.svg")

# %% Stacked plot by percentage


def adjustFigAspect(fig, aspect:float=1):
    """
    Adjust subplot parameters to achieve the correct aspect ratio.

    This function modifies the figure layout to ensure proper proportions
    for publication-quality plots.

    Args:
        fig (matplotlib.figure.Figure): Figure object to adjust
        aspect (float): Desired aspect ratio (default: 1)
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
    """
    Generate a stacked bar chart showing assay results by pathogen genus.

    Creates a visualization of epitope binding assay results grouped by
    pathogen genus, with bars showing the distribution of qualitative
    measures (positive, negative, etc.) for each genus.

    Args:
        data (pd.DataFrame): Merged dataset with IEDB and taxonomy information
        division (str): Taxonomic division to filter by (e.g., "Bacteria", "Viruses")
        output_file (str): Path for the output SVG file

    Note:
        The function filters data by the specified division and creates
        a stacked percentage bar chart showing assay result distributions.
    """

    grouped = (
        data[data["Division"] == division]
        .groupby(["genus", "qualitative_measure"])
        .count()
    )
    grouped = grouped["assay_iris"].unstack()
    n_samples = grouped.sum(axis=1)
    grouped.index = pd.Index([
        "{} (n={})".format(y.replace("_", " ").title(), int(val))
        for y, val in n_samples.to_dict().items()
    ])
    grouped = grouped[grouped.sum(axis=1) > 10]

    # sort by the sum of the columns
    grouped = grouped.loc[grouped.sum(
        axis=1).sort_values(ascending=False).index]

    # normalize rows
    grouped = grouped.div(grouped.sum(axis=1), axis=0)

    colors = {
        "Negative": "#8b0000",
        "Positive-Low": "#f6655f",
        "Positive": "#c1f57d",
        "Positive-Intermediate": "#97e692",
        "Positive-High": "#069d59"
    }
    # Filter qualitative values
    colors = {key: val for key, val in colors.items() if key in grouped.columns}

    # Sort columns
    qual_measurements = list(colors.keys())
    grouped = grouped.loc[:, qual_measurements]

    # Invert negative values
    all_negatives = ["Negative", "Positive-Low"]
    negative_values = [val for val in grouped.columns if val in all_negatives]
    grouped[negative_values] = grouped[negative_values] * -1

    fig, ax = plt.subplots()
    adjustFigAspect(fig, aspect=0.7)
    grouped.plot(kind="barh", stacked=True, width=0.7, color=colors, ax=ax)
    plt.legend(loc="upper left", bbox_to_anchor=(1, 1))
    # set the x-axis limits to -1, 1
    ax.set_xlim(-1, 1)
    ax.autoscale(enable=True, axis="y", tight=True)

    plt.savefig(output_file, bbox_inches="tight")
