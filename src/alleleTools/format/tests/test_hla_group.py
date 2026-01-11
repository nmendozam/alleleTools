import pytest

from ..hla_group import GrouperHLA

def test_3f_allele_grouping():
    grouper = GrouperHLA("g-group")

    group = grouper.lookup("DQB1", "02:01:01")
    assert group == "DQB1*02:01:01G"

    group = grouper.lookup("DQB1", "02:02:01")
    assert group == "DQB1*02:01:01G"

def test_2f_allele_grouping():
    grouper = GrouperHLA("g-group")

    group = grouper.lookup("DQB1", "02:01")
    assert group == "DQB1*02:01G"

    group = grouper.lookup("DQB1", "02:02")
    assert group == "DQB1*02:01G"