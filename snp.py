from __future__ import print_function

from decimal import Decimal


class SnpFeature:
    GENE_LENGTH = 2

    def __init__(self, snp):
        features = SnpFeature.extract_feature(snp)
        self.data = features
        self.large_count = features.count(1)
        self.small_count = features.count(-1)
        self.medium_count = len(features) - self.large_count - self.small_count

    @staticmethod
    def extract_feature(snp):
        snp00 = snp[0][0]
        large_item = small_item = snp00
        for item in snp[0]:
            if item > large_item:
                large_item = item
            if item < small_item:
                small_item = item

        return [SnpFeature.__extract_feature_from_symbol(snp[0][i] + snp[1][i], large_item, small_item)
                for i in range(0, len(snp[0]))]

    @staticmethod
    def __extract_feature_from_symbol(symbol, large_item, small_item):
        if symbol == large_item * 2:
            return 1
        elif symbol == small_item * 2:
            return -1
        else:
            return 0

    def p_A(self):
        large_count = self.large_count * self.GENE_LENGTH + self.medium_count

        return Decimal(large_count) / (len(self.data) * self.GENE_LENGTH)

    def p_AA_Aa_aa(self):
        return tuple(map(lambda x: Decimal(x) / len(self.data),
                         (self.large_count, self.medium_count, self.small_count)))
