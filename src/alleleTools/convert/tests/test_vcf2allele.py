from ..vcf2allele import VCF
from ...utils.assets import get_asset_path


def test_read_vcf():
    vcf_file = get_asset_path("kir_example.vcf")
    vcf = VCF(vcf_file)

    expected_columns = [
        "CHROM",
        "POS",
        "REF",
        "ALT",
        "QUAL",
        "FILTER",
        "INFO",
        "FORMAT",
    ]
    first_columns = vcf.dataframe.columns.to_list()[: len(expected_columns)]
    assert first_columns == expected_columns

    assert vcf.dataframe.index.name == "ID"
