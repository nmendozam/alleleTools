"""
PheWAS Visualization Module.

This module provides functionality to create PheWAS (Phenome-Wide Association Study)
Manhattan plots for HLA alleles. It processes data from the PheWAS catalog and
generates publication-quality visualizations showing disease associations across
different phenotype categories.

Data source: https://phewascatalog.org/phewas/#hla
The input file should be downloaded from the PheWAS catalog.

Author: Nicolás Mendoza Mejía (2025)
"""

from typing import Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ..argtypes import file_path, output_path


def setup_parser(subparsers):
    """
    Set up the argument parser for the pleiotropy command.

    Args:
        subparsers: The subparsers object to add this command to.

    Returns:
        argparse.ArgumentParser: The configured parser for pleiotropy visualization.
    """
    parser = subparsers.add_parser(
        name="pleiotropy",
        help="Crete pleiotropy report",
        description="Create peliotropy report (Manhattan graph of phewas)",
        epilog="Author: Nicolás Mendoza Mejía (2025)",
    )
    parser.add_argument(
        "input",
        type=file_path,
        help="path to the CSV downlowaded from https://phewascatalog.org/phewas/#hla",
    )
    parser.add_argument(
        "--output",
        type=output_path,
        help="name of the output image file generated",
        default="pleiotropy.svg"
    )
    parser.add_argument(
        "--allele_name",
        type=str,
        help="Name of the allele to include in the title",
        default=""
    )

    parser.set_defaults(func=call_function)

    return parser


def call_function(args):
    """
    Execute PheWAS plot generation.

    Args:
        args: Parsed command line arguments containing:
            - input: Path to PheWAS catalog CSV file
            - output: Path to output image file
            - allele_name: Name of allele for plot title
    """
    create_phewas_plot(args.input, args.output, allele_name=args.allele_name)


def create_phewas_plot(input_file: str,
                       output_file: str,
                       figsize: Tuple[int, ...] = (15, 6),
                       significance_threshold: float = 10**-5,
                       allele_name: str = '',
                       ):
    """
    Create a PheWAS Manhattan plot from PheWAS catalog data.

    Generates a Manhattan plot showing disease associations for an HLA allele
    across different phenotype categories. Points are colored by category and
    sized by significance level.

    Args:
        input_file (str): Path to CSV file from PheWAS catalog
        output_file (str): Path for output image file
        figsize (Tuple[int, ...]): Figure size in inches (default: (15, 6))
        significance_threshold (float): P-value threshold for significance line (default: 1e-5)
        allele_name (str): Allele name to include in plot title (default: '')

    Note:
        Input CSV should contain columns: 'phecode', 'phenotype', 'p', 'category'
        as provided by the PheWAS catalog (https://phewascatalog.org/phewas/#hla)
    """

    try:
        df = pd.read_csv(input_file)
    except Exception as e:
        raise Exception(
            f"An error occurred while opening the file {input_file}: {e}")

    # filter by allele name
    if allele_name:
        df = df[df['HLA'].str.contains(allele_name, case=False, na=False)]
        if df.empty:
            raise ValueError(
                f"No data found for allele name '{allele_name}' in the input file. Check that the argument --allele_name matches with the provided input file.")

    df.sort_values(by='Category', inplace=True)
    df.index = range(len(df.index))
    df["position"] = df.index

    # Calculate -log10(P-value)
    df['-log10(P)'] = -np.log10(df['P-Value'])

    # Create a dictionary to map categories to colors
    group_cat = df.reset_index().groupby('Category')
    categories = group_cat['Category'].first().unique()
    # Use a colormap for distinct colors
    colors = plt.cm.get_cmap('tab20', len(categories))
    category_colors = {cat: colors(i) for i, cat in enumerate(categories)}

    # Create the Manhattan plot
    plt.figure(figsize=figsize)

    # Plot each point, coloring by category
    for category, color in category_colors.items():
        subset = df[df['Category'] == category]
        # Plot the points for this category with alpha 0.8 and a contour line
        plt.scatter(subset["position"], subset['-log10(P)'],
                    color=color, label=category, s=20, alpha=0.5)

    # Set x-ticks to be the Category names
    # place them at the middle of the category
    beginning = group_cat['position'].min()
    sizes = group_cat['position'].size()
    middle = ((sizes / 2) + beginning).sort_values()
    plt.xticks(middle.to_list(), middle.index, rotation=45, ha='right')

    # Set x and y limits
    plt.xlim(0, len(df.index))  # Adjust x-axis limits based on your data
    # Add some space above the highest point
    plt.ylim(0, df['-log10(P)'].max() + 1)

    # Customize the plot
    plt.ylabel('-log10(P-value)')
    plt.title(f"PheWAS Manhattan Plot {allele_name}")
    # Rotate x-axis labels for readability
    plt.xticks(rotation=45, fontsize=10)

    # color of xticks by category
    for tick in plt.gca().get_xticklabels():
        tick.set_color(category_colors[tick.get_text()])

    # Add labels for significant points (optional)
    # Adjust significance threshold as needed
    significant_points = df[df['P-Value'] < 10**-5]
    for index, row in significant_points.iterrows():
        plt.annotate(row['Phenotype'], (row['position'], row['-log10(P)']),
                     textcoords="offset points", xytext=(0, 5), ha='center', fontsize=8)

    # Add a horizontal line for significance threshold (optional)
    if significance_threshold >= 0:
        plt.axhline(y=5, color='r', linestyle='--',
                    label='Significance (%f)' % 10**-5)

    # Add legend
    plt.legend(loc='upper left', bbox_to_anchor=(1, 1))

    # Show the plot
    plt.tight_layout()  # Adjust layout to prevent labels from overlapping
    # plt.show()

    # save as svg
    plt.savefig(output_file, format='svg', bbox_inches='tight')
