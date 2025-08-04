import pandas as pd

from .ukb2allele import (_add_missing_alleles, _convert_ukb_to_allele,
                         _format_allele_names)


def test_convert_ukb_to_allele():
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

    case_control = _convert_ukb_to_allele(input, phenotype)

    expected = pd.DataFrame(
        {
            "eid": [1, 2, 3],
            "Pheno": [1, 2, 0],
            "A": ["A*01:01", "A*01:02", "A*01:01"],
            "A_2": ["A*01:02", "NA", "NA"],
            "B": ["B*02:01", "B*02:01", "NA"],
            "B_2": ["NA", "NA", "NA"],
        }
    )

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

    ret = _add_missing_alleles(row)
    assert ret == "A_101,A_NA,B_201,B_NA"


def test_fill_allele_pairs_unsorted():
    row = "B_201,A_101,A_102"

    ret = _add_missing_alleles(row)
    assert ret == "A_101,A_102,B_201,B_NA"


def test_allele_formatting():
    df = pd.DataFrame(
        {
            "A": ["A_101", "A_NA"],
            "A_2": ["A_102", "A_NA"],
        }
    )

    formatted = _format_allele_names(df)

    expected = {
        "A": ["A*01:01", "NA"],
        "A_2": ["A*01:02", "NA"],
    }

    pd.testing.assert_frame_equal(formatted, pd.DataFrame(expected))
