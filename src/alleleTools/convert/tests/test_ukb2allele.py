import pandas as pd

from ..ukb2allele import (__format_allele_names, _convert_ukb_to_allele,
                          _na_missing_alleles)


def test_convert_ukb_to_allele():
    phenotype = {
        "eid": [1, 2, 3],
        "Pheno": [1, 2, 0],
    }

    input = {
        "eid": [1, 2, 3],
        "A_101": [1, 0, 1],
        "B_201": [1, 1, 0],
        "C_101": [1, 0, 1],
        "DRB5_101": [1, 0, 1],
        "DRB4_101": [1, 0, 1],
        "DRB3_101": [1, 0, 1],
        "DRB1_101": [1, 0, 1],
        "DQB1_101": [1, 0, 1],
        "DQA1_101": [1, 0, 1],
        "DPB1_101": [1, 0, 1],
        "DPA1_101": [1, 0, 1],
    }

    input = pd.DataFrame(input)
    input.set_index("eid", inplace=True)

    phenotype = pd.DataFrame(phenotype)

    case_control = _convert_ukb_to_allele(input, phenotype)

    expected = pd.DataFrame(
        {
            "eid": [1, 2, 3],
            "Pheno": [1, 2, 0],
            "A": ["A*01:01", "NA", "A*01:01"],
            "A_2": ["NA", "NA", "NA"],
            "B": ["B*02:01", "B*02:01", "NA"],
            "B_2": ["NA", "NA", "NA"],
            "C": ["C*01:01", "NA", "C*01:01"],
            "C_2": ["NA", "NA", "NA"],
            "DRB5": ["DRB5*01:01", "NA", "DRB5*01:01"],
            "DRB5_2": ["NA", "NA", "NA"],
            "DRB4": ["DRB4*01:01", "NA", "DRB4*01:01"],
            "DRB4_2": ["NA", "NA", "NA"],
            "DRB3": ["DRB3*01:01", "NA", "DRB3*01:01"],
            "DRB3_2": ["NA", "NA", "NA"],
            "DRB1": ["DRB1*01:01", "NA", "DRB1*01:01"],
            "DRB1_2": ["NA", "NA", "NA"],
            "DQB1": ["DQB1*01:01", "NA", "DQB1*01:01"],
            "DQB1_2": ["NA", "NA", "NA"],
            "DQA1": ["DQA1*01:01", "NA", "DQA1*01:01"],
            "DQA1_2": ["NA", "NA", "NA"],
            "DPB1": ["DPB1*01:01", "NA", "DPB1*01:01"],
            "DPB1_2": ["NA", "NA", "NA"],
            "DPA1": ["DPA1*01:01", "NA", "DPA1*01:01"],
            "DPA1_2": ["NA", "NA", "NA"],
        }
    )
    print(case_control)
    print(expected)

    expected = expected.reindex(sorted(expected.columns), axis=1)
    case_control = case_control.reindex(sorted(case_control.columns), axis=1)

    pd.testing.assert_frame_equal(case_control, expected)


def test_remove_pheno_zero():
    phenotype = {
        "eid": [1, 2, 3],
        "Pheno": [1, 2, 0],
    }

    input = {
        "eid": [1, 2, 3],
        "A_101": [1, 0, 1],
        "A_102": [1, 1, 0],
        "B_201": [1, 1, 0],
    }

    input = pd.DataFrame(input)
    input.set_index("eid", inplace=True)

    phenotype = pd.DataFrame(phenotype)

    case_control = _convert_ukb_to_allele(input, phenotype, rm_phe_zero=True)

    expected = pd.DataFrame(
        {
            "eid": [1, 2],
            "Pheno": [1, 2],
            "A": ["A*01:01", "A*01:02"],
            "A_2": ["A*01:02", "NA"],
            "B": ["B*02:01", "B*02:01"],
            "B_2": ["NA", "NA"],
        }
    )

    pd.testing.assert_frame_equal(case_control, expected)


def test_fill_allele_pairs():
    row = "A_101,B_201"

    ret = _na_missing_alleles(row)
    assert ret == "A_101,A_NA,B_201,B_NA"


def test_fill_allele_pairs_unsorted():
    row = "B_201,A_101,A_102"

    ret = _na_missing_alleles(row)
    assert ret == "A_101,A_102,B_201,B_NA"


def test_allele_formatting():
    df = pd.DataFrame(
        {
            "A": ["A_101", "A_NA"],
            "A_2": ["A_102", "A_NA"],
        }
    )

    formatted = __format_allele_names(df)

    expected = {
        "A": ["A*01:01", "NA"],
        "A_2": ["A*01:02", "NA"],
    }

    pd.testing.assert_frame_equal(formatted, pd.DataFrame(expected))
