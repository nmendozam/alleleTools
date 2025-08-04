from typing import Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from ..argtypes import file_path, output_path


def setup_parser(subparsers):
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
    create_phewas_plot(args.input, args.output, allele_name=args.allele_name)

def create_phewas_plot(input_file: str, 
                       output_file: str,
                       figsize: Tuple[int, ...] = (15, 6),
                       significance_threshold: float = 10**-5,
                       allele_name: str = '',
                       ):
    """
    Creates a PheWAS Manhattan plot from a CSV file.
    """

    try:
        df = pd.read_csv(input_file)
    except Exception as e:
        raise Exception(f"An error occurred while opening the file {input_file}: {e}")
    
    # filter by allele name
    if allele_name:
        df = df[df['HLA'].str.contains(allele_name, case=False, na=False)]
        if df.empty:
            raise ValueError(f"No data found for allele name '{allele_name}' in the input file. Check that the argument --allele_name matches with the provided input file.")

    df.sort_values(by='Category', inplace=True)
    df.index = range(len(df.index))
    df["position"] = df.index

    # Calculate -log10(P-value)
    df['-log10(P)'] = -np.log10(df['P-Value'])

    # Create a dictionary to map categories to colors
    group_cat= df.reset_index().groupby('Category')
    categories = group_cat['Category'].first().unique()
    colors = plt.cm.get_cmap('tab20', len(categories))  # Use a colormap for distinct colors
    category_colors = {cat: colors(i) for i, cat in enumerate(categories)}

    # Create the Manhattan plot
    plt.figure(figsize=figsize)

    # Plot each point, coloring by category
    for category, color in category_colors.items():
        subset = df[df['Category'] == category]
        # Plot the points for this category with alpha 0.8 and a contour line
        plt.scatter(subset["position"], subset['-log10(P)'], color=color, label=category, s=20, alpha=0.5)
    
    # Set x-ticks to be the Category names
    # place them at the middle of the category
    beginning = group_cat['position'].min()
    sizes = group_cat['position'].size()
    middle = ((sizes / 2) + beginning).sort_values()
    plt.xticks(middle.to_list(), middle.index, rotation=45, ha='right')

    # Set x and y limits
    plt.xlim(0, len(df.index))  # Adjust x-axis limits based on your data
    plt.ylim(0, df['-log10(P)'].max() + 1)  # Add some space above the highest point


    # Customize the plot
    plt.ylabel('-log10(P-value)')
    plt.title(f"PheWAS Manhattan Plot {allele_name}")
    plt.xticks(rotation=45, fontsize=10)  # Rotate x-axis labels for readability

    # color of xticks by category
    for tick in plt.gca().get_xticklabels():
        tick.set_color(category_colors[tick.get_text()])

    # Add labels for significant points (optional)
    significant_points = df[df['P-Value'] < 10**-5]  # Adjust significance threshold as needed
    for index, row in significant_points.iterrows():
        plt.annotate(row['Phenotype'], (row['position'], row['-log10(P)']), 
                        textcoords="offset points", xytext=(0, 5), ha='center', fontsize=8)

    # Add a horizontal line for significance threshold (optional)
    if significance_threshold >= 0:
        plt.axhline(y=5, color='r', linestyle='--', label='Significance (%f)' % 10**-5)

    # Add legend
    plt.legend(loc='upper left', bbox_to_anchor=(1, 1))


    # Show the plot
    plt.tight_layout()  # Adjust layout to prevent labels from overlapping
    # plt.show()

    # save as svg
    plt.savefig(output_file, format='svg', bbox_inches='tight')

