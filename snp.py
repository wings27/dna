from __future__ import print_function

from decimal import Decimal
import numpy


class FileHelper:
    @staticmethod
    def load_feature_group(file_name):
        return numpy.loadtxt(file_name, int, delimiter=',')


class SnpFeature:
    GENE_LENGTH = 2

    def __init__(self, features):
        self.data = features
        self.__largeCount = features.count(1)
        self.__smallCount = features.count(-1)
        self.__mediumCount = len(features) - self.__largeCount - self.__smallCount

    def __getitem__(self, item):
        return self.data[item]

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
        large_count = self.__largeCount * self.GENE_LENGTH + self.__mediumCount

        return Decimal(large_count) / (len(self.data) * self.GENE_LENGTH)

    def p_AA_Aa_aa(self):
        return tuple(map(lambda x: Decimal(x) / len(self.data),
                         (self.__largeCount, self.__mediumCount, self.__smallCount)))


class CLDCalculation:
    def __init__(self, snp_group_1, snp_group_2):
        self.snp_group_1 = snp_group_1
        self.snp_group_2 = snp_group_2
        self.n = len(snp_group_1[0])

    def cld_1(self):
        cld1_max = -100000
        for line1 in self.snp_group_1:
            for line2 in self.snp_group_2:
                f1 = SnpFeature(list(line1))
                f2 = SnpFeature(list(line2))
                n_AB_value = self.__n_AB(f1, f2)
                sig_AB = n_AB_value / self.n - 2 * f1.p_A() * f2.p_A()
                down_left = f1.p_A() * (1 - f1.p_A()) + f1.p_AA_Aa_aa()[0] - f1.p_A() * f1.p_A()
                down_right = f2.p_A() * (1 - f2.p_A()) + f2.p_AA_Aa_aa()[0] - f2.p_A() * f2.p_A()

                cld1 = self.n * sig_AB ** 2 / (down_left * down_right)
                if cld1 > cld1_max:
                    cld1_max = cld1
                print(cld1)
        return cld1_max

    def __n_AB(self, snp_feature1, snp_feature2):
        n_count = {(1, 1): 0, (1, 0): 0, (0, 1): 0, (0, 0): 0}
        for i in range(0, self.n):
            key = (snp_feature1[i], snp_feature2[i])
            if key in n_count:
                n_count[key] += 1
        return 2 * n_count[(1, 1)] + n_count[(1, 0)] + n_count[(0, 1)] + Decimal(0.5) * n_count[(0, 0)]
