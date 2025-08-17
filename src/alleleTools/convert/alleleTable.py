import pandas as pd


class AlleleTable:
    def __init__(self):
        self.alleles = pd.DataFrame()
        self.phenotype = pd.Series()
        self.covariates = pd.DataFrame()
