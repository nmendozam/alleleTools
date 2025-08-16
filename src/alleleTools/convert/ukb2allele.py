# %%
import sys

import pandas as pd

from ..argtypes import csv_file, output_path


def setup_parser(subparsers):
    parser = subparsers.add_parser(
        name="ukb2allele",
        description="Convert UK Biobank HLA data to allele table format",
        epilog="Author: Nicolás Mendoza Mejía (2025)",
    )
    ## Input/output arguments
    parser.add_argument(
        "input",
        type=csv_file,
        help="Input csv file with the UK Biobank imputed HLA data",
    )
    parser.add_argument(
        "--phenotype",
        type=csv_file,
        help="ssv file with 6 columns: eid, fid, ... , Sex, Pheno. No headers and space separated. The column Pheno (last column) will be included as phenotype in the output file",
        required=True,
    )
    parser.add_argument(
        "--output",
        type=output_path,
        help="name of the output file",
        default="output.alt",
    )
    ## Additional arguments
    parser.add_argument(
        "--remove_pheno_zero",
        action="store_true",
        help="Remove individuals with phenotype 0 from the output",
        default=False,
    )

    parser.set_defaults(func=call_function)

    return parser


def call_function(args):
    input = pd.read_csv(args.input, sep=",", header=0, index_col=0)

    phenotype = pd.read_csv(args.phenotype, sep=" ", header=None)
    phenotype.columns = ["eid", "FID", "_", "__", "Sex", "Pheno"]

    output = _convert_ukb_to_allele(input, phenotype, args.remove_pheno_zero)

    # Save the result to a file
    output.to_csv(args.output, sep="\t", index=False, header=False)

    print(output.info())
    print(output["Pheno"].value_counts())
    print(f"Output saved to {args.output}")


def __get_list_of_genes(df: pd.DataFrame) -> list:
    return df.columns.map(lambda x: x.split("_")[0]).unique().tolist()


def __filter_alleles(
    df: pd.DataFrame, min_probability: float = 0.7, homo_thr: float = 1.4
) -> pd.DataFrame:
    """
    Filters alleles based on their probability. The UK Biobank recommends filtering out alleles with less than 70% probability. Aleles with a probability greater than 1.4 are considered homozygous and duplicated in the output.
    - df = dataframe with the allele table
    - min_probability = minimum probability threshold
    - homo_thr = homozygous threshold
    """
    # Get only values with non-zero probability
    df_filtered = df[df["Presence"] > 0]
    # Filter out alleles with less than 70% probability as recommended by UKB
    df_filtered.loc[
        df_filtered["Presence"] <= min_probability, "Allele"
    ] = df_filtered.loc[df_filtered["Presence"] <= min_probability, "Allele"].apply(
        lambda x: x.split("_")[0] + "_NA"
    )
    # Duplicate homozygous alleles
    df_filtered.loc[df_filtered["Presence"] > homo_thr, "Allele"] = (
        df_filtered.loc[df_filtered["Presence"] > homo_thr, "Allele"]
        + ","
        + df_filtered.loc[df_filtered["Presence"] > homo_thr, "Allele"]
    )

    return df_filtered


def __assign_genes(str):
    alleles = dict()
    for allele in str.split(","):
        locus = allele.split("_")[0]
        if locus in alleles.keys():
            locus += "_2"
        alleles[locus] = allele

    return pd.Series(alleles.values(), index=alleles.keys())


def _convert_ukb_to_allele(
    input: pd.DataFrame, phenotype: pd.DataFrame, rm_phe_zero: bool = False
) -> pd.DataFrame:
    df_melted = input.reset_index().melt(
        id_vars="eid", var_name="Allele", value_name="Presence"
    )

    # Filter rows where the allele is present
    df_filtered = __filter_alleles(df_melted)

    # Joint allele names from each individual
    joint_alleles = df_filtered.groupby("eid")["Allele"].apply(",".join).reset_index()

    # Fills in the gaps of missing alleles
    joint_alleles["Allele"] = joint_alleles["Allele"].apply(_na_missing_alleles)

    # Reshape the table
    expanded = joint_alleles.Allele.apply(__assign_genes)

    # Concatenate individual ID with alleles
    df_alleles = pd.concat([joint_alleles["eid"], expanded], axis=1)

    # Properly format allele names
    df_formatted_alleles = _format_allele_names(df_alleles)

    # Add the case-control values
    # Get the Pheno values from phenotype and assign it to a new column in df_formatted_alleles
    df_case_control = df_formatted_alleles.merge(
        phenotype[["eid", "Pheno"]], on="eid", how="inner"
    )

    # place it on the second column
    df_case_control = df_case_control[
        ["eid", "Pheno"]
        + [col for col in df_case_control.columns if col != "Pheno" and col != "eid"]
    ]
    df_case_control.head()

    # remove Pheno with 0
    if rm_phe_zero:
        df_case_control = df_case_control[df_case_control["Pheno"] != 0]

    return df_case_control


def _format_allele_names(expanded: pd.DataFrame) -> pd.DataFrame:
    """
    Reformats allele names to the standard format:
    - A_101 becomes A*01:01
    - Fills missing alleles with NA
    """
    # Rename columns to the desired format
    df = expanded.replace(to_replace=r".*_NA", value="NA", regex=True)
    df = df.replace(to_replace=r".*9901", value="NA", regex=True)
    df = df.fillna("NA")
    # Replace all _ with *
    df = df.replace(to_replace=r"_", value="*", regex=True)
    # Add : before the last 2 digits
    df = df.replace(to_replace=r"([0-9]{2})\b", value=":\\1", regex=True)
    # if there is less than 4 digits after the *, add a 0 at the beginning
    df = df.replace(to_replace=r"\*([0-9]{1}):", value="*0\\1:", regex=True)
    # Remove any spaces
    df = df.replace(to_replace=r"\s", value="", regex=True)
    return df


def _na_missing_alleles(row):
    """
    This expects a list of alleles as a string separated by commas.
    It returns a list of alleles with NA added for missing pairs.
    For example, if the input is "A_101,B_201", the output will
    be "A_101,A_NA,B_201,B_NA".
    """
    alleles = row.split(",")

    # Sort alleles
    alleles.sort()

    checked_alleles = []
    # count how many duplicated genes there are
    current = alleles[0].split("_")[0]
    num_alleles = 0
    for allele in alleles:
        next = allele.split("_")[0]

        if next == current:
            num_alleles += 1
        else:  # Change of new gene
            if num_alleles < 2:
                checked_alleles.append(current + "_NA")
            current = next
            num_alleles = 1

        checked_alleles.append(allele)

    # If the last gene has less than 2 alleles, add NA
    if len(checked_alleles) % 2 != 0:
        checked_alleles.append(current + "_NA")

    return ",".join(checked_alleles)
