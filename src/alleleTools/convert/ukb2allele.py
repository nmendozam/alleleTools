#%%
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
    input = pd.read_csv(args.input, sep=',', header=0, index_col=0)

    phenotype = pd.read_csv(args.phenotype, sep=' ', header=None)
    phenotype.columns = ["eid", "FID", "_", "__", "Sex", "Pheno"]

    output = _convert_ukb_to_allele(input, phenotype, args.remove_pheno_zero)

    # Save the result to a file
    output.to_csv(args.output, sep='\t', index=False, header=False)

    print(output.info())
    print(output["Pheno"].value_counts())
    print(f"Output saved to {args.output}")

def _convert_ukb_to_allele(input: pd.DataFrame, phenotype: pd.DataFrame, rm_phe_zero: bool = False) -> pd.DataFrame:
    #%%
    df_melted = input.reset_index().melt(id_vars="eid", var_name="Allele", value_name="Presence")
    df_melted.head()

    #%% Filter rows where the allele is present

    # Get only values with non-zero probability
    df_filtered = df_melted[df_melted["Presence"] > 0]
    # Filter out alleles with less than 70% probability as recommended by UKB
    df_filtered.loc[df_filtered["Presence"] <= 0.7, "Allele"] = df_filtered.loc[df_filtered["Presence"] <= 0.7, "Allele"].apply(lambda x: x.split("_")[0] + "_NA")
    # Duplicate homozigous alleles
    df_filtered.loc[df_filtered["Presence"] > 1.4, "Allele"] = df_filtered.loc[df_filtered["Presence"] > 1.4, "Allele"] + "," + df_filtered.loc[df_filtered["Presence"] > 1.4, "Allele"]


    #%% Joint allele names from each individual
    joint_alleles = df_filtered.groupby("eid")["Allele"].apply(",".join).reset_index()


    # For the alleles with less than 22 items, add NA to make it 22
    joint_alleles["Allele"] = joint_alleles["Allele"].apply(_add_missing_alleles)

    #%%
    # Expand alleles into individual columns and name them
    expanded = joint_alleles["Allele"].str.split(',', expand=True)
    # get number of columns
    num_cols = expanded.shape[1]
    expanded.columns = ["A", "A_2", "B", "B_2", "C", "C_2", "DRB5", "DRB5_2", "DRB4", "DRB4_2", "DRB3", "DRB3_2", "DRB1", "DRB1_2", "DQB1", "DQB1_2", "DQA1", "DQA1_2", "DPB1", "DPB1_2", "DPA1", "DPA1_2"][0:num_cols]
    # Concatenate individual ID with alleles
    df_alleles = pd.concat([joint_alleles["eid"], expanded], axis=1)


    #%% Properly format allele names
    df_formatted_alleles = _format_allele_names(df_alleles)

    #%% Add the case-control values
    # Get the Pheno values from phenotype and assign it to a new column in df_formatted_alleles
    df_case_control = df_formatted_alleles.merge(phenotype[["eid", "Pheno"]], on="eid", how="inner")

    # place it on the second column
    df_case_control = df_case_control[["eid", "Pheno"] + [col for col in df_case_control.columns if col != "Pheno" and col != "eid"]]
    df_case_control.head()

    # remove Pheno with 0
    if rm_phe_zero:
        df_case_control = df_case_control[df_case_control["Pheno"] != 0]

    return df_case_control

def _format_allele_names(expanded: pd.DataFrame) -> pd.DataFrame:
    # Rename columns to the desired format
    df = expanded.replace(to_replace=r'.*_NA', value='NA', regex=True)
    df = df.replace(to_replace=r'.*9901', value='NA', regex=True)
    df = df.fillna('NA')
    # Replace all _ with *
    df = df.replace(to_replace=r'_', value='*', regex=True)
    # Add : before the last 2 digits
    df = df.replace(to_replace=r'([0-9]{2})\b', value=':\\1', regex=True)
    # if there is less than 4 digits after the *, add a 0 at the beginning
    df = df.replace(to_replace=r'\*([0-9]{1}):', value='*0\\1:', regex=True)
    # Remove any spaces
    df = df.replace(to_replace=r'\s', value='', regex=True)
    return df

def _add_missing_alleles(row):
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
        else: # Change of new gene
            if num_alleles < 2:
                checked_alleles.append(current + "_NA")
            current = next
            num_alleles = 1

        checked_alleles.append(allele)

    # If the last gene has less than 2 alleles, add NA
    if len(checked_alleles)%2 != 0:
        checked_alleles.append(current + "_NA")

    return ",".join(checked_alleles)
