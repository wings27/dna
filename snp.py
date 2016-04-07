from __future__ import print_function

from decimal import Decimal

import numpy


class FileHelper:
    @staticmethod
    def load_feature_group(file_name):
        loaded_result = numpy.loadtxt(file_name, int, delimiter=',')
        if len(loaded_result.shape) == 1:
            loaded_result = loaded_result.reshape((1, -1))
        return loaded_result

    @staticmethod
    def save_array(file_name, array, fmt):
        with open(file_name + '.out', 'ab') as f:
            numpy.savetxt(f, array, fmt, ',')


class SnpFeature:
    GENE_LENGTH = 2
    EMPTY = 2

    def __init__(self, features):
        self.data = features
        self.__largeCount = features.count(1)
        self.__smallCount = features.count(-1)
        self.__mediumCount = len(features) - self.__largeCount - self.__smallCount

    def __getitem__(self, item):
        return self.data[item]

    @staticmethod
    def extract_feature(snp):
        chars = set(snp[0]).union(set(snp[1])).difference(['0'])
        small_item = sorted(chars)[0]
        large_item = sorted(chars)[-1]

        return [SnpFeature.__extract_feature_from_symbol(snp[0][i] + snp[1][i], large_item, small_item)
                for i in range(len(snp[0]))]

    @staticmethod
    def __extract_feature_from_symbol(symbol, large_item, small_item):
        if symbol == large_item * 2:
            return 1
        elif symbol == small_item * 2:
            return -1
        elif symbol == '00':
            return SnpFeature.EMPTY
        else:
            return 0

    def p_A(self):
        large_count = self.__largeCount * self.GENE_LENGTH + self.__mediumCount

        return Decimal(large_count) / ((len(self.data) - self.data.count(SnpFeature.EMPTY)) * self.GENE_LENGTH)

    def p_AA_Aa_aa(self):
        return tuple(map(lambda x: Decimal(x) / (len(self.data) - self.data.count(SnpFeature.EMPTY)),
                         (self.__largeCount, self.__mediumCount, self.__smallCount)))


class CLDCalculation:
    F_AABB = (1, 1)
    F_AABb = (1, 0)
    F_AAbb = (1, -1)
    F_AaBB = (0, 1)
    F_AaBb = (0, 0)
    F_Aabb = (0, -1)
    F_aaBB = (-1, 1)
    F_aaBb = (-1, 0)
    F_aabb = (-1, -1)

    def __init__(self, snp_group_1, snp_group_2):
        self.snp_group_1 = snp_group_1
        self.snp_group_2 = snp_group_2
        self.n = len(snp_group_1[0])

    def temp_AB(self):
        temp_max = -100000
        for line1 in self.snp_group_1:
            for line2 in self.snp_group_2:
                f1 = SnpFeature(list(line1))
                f2 = SnpFeature(list(line2))
                n_AB_value = self.__n_AB(f1, f2)
                temp = n_AB_value / self.n - 2 * f1.p_A() * f2.p_A()

                if temp > temp_max:
                    temp_max = temp
        return temp_max

    def temp_Ab(self):
        temp_max = -100000
        for line1 in self.snp_group_1:
            for line2 in self.snp_group_2:
                f1 = SnpFeature(list(line1))
                f2 = SnpFeature(list(line2))
                n_Ab_value = self.__n_Ab(f1, f2)
                temp = n_Ab_value / self.n - 2 * f1.p_A() * (1 - f2.p_A())

                if temp > temp_max:
                    temp_max = temp
        return temp_max

    def temp_aB(self):
        temp_max = -100000
        for line1 in self.snp_group_1:
            for line2 in self.snp_group_2:
                f1 = SnpFeature(list(line1))
                f2 = SnpFeature(list(line2))
                n_aB_value = self.__n_aB(f1, f2)
                temp = n_aB_value / self.n - 2 * (1 - f1.p_A()) * f2.p_A()

                if temp > temp_max:
                    temp_max = temp
        return temp_max

    def temp_ab(self):
        temp_max = -100000
        for line1 in self.snp_group_1:
            for line2 in self.snp_group_2:
                f1 = SnpFeature(list(line1))
                f2 = SnpFeature(list(line2))
                n_ab_value = self.__n_ab(f1, f2)
                temp = n_ab_value / self.n - 2 * (1 - f1.p_A()) * (1 - f2.p_A())

                if temp > temp_max:
                    temp_max = temp
        return temp_max

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
        return cld1_max

    def cld_2(self):
        cld2_max = -100000
        for line1 in self.snp_group_1:
            for line2 in self.snp_group_2:
                f1 = SnpFeature(list(line1))
                f2 = SnpFeature(list(line2))
                n_c = self.__n_count(f1, f2)
                p_AA_x_p_BB = Decimal(f1.p_AA_Aa_aa()[0]) * f2.p_AA_Aa_aa()[0]
                adder_1 = (Decimal(n_c.get(CLDCalculation.F_AABB, 0)) / self.n - p_AA_x_p_BB) ** 2 / p_AA_x_p_BB

                p_AA_x_p_Bb = Decimal(f1.p_AA_Aa_aa()[0]) * f2.p_AA_Aa_aa()[1]
                adder_2 = (Decimal(n_c.get(CLDCalculation.F_AABb, 0)) / self.n - p_AA_x_p_Bb) ** 2 / p_AA_x_p_Bb

                p_AA_x_p_bb = Decimal(f1.p_AA_Aa_aa()[0]) * f2.p_AA_Aa_aa()[2]
                adder_3 = (Decimal(n_c.get(CLDCalculation.F_AAbb, 0)) / self.n - p_AA_x_p_bb) ** 2 / p_AA_x_p_bb

                p_Aa_x_p_BB = Decimal(f1.p_AA_Aa_aa()[1]) * f2.p_AA_Aa_aa()[0]
                adder_4 = (Decimal(n_c.get(CLDCalculation.F_AaBB, 0)) / self.n - p_Aa_x_p_BB) ** 2 / p_Aa_x_p_BB

                p_Aa_x_p_Bb = Decimal(f1.p_AA_Aa_aa()[1]) * f2.p_AA_Aa_aa()[1]
                adder_5 = (Decimal(n_c.get(CLDCalculation.F_AaBb, 0)) / self.n - p_Aa_x_p_Bb) ** 2 / p_Aa_x_p_Bb

                p_Aa_x_p_bb = Decimal(f1.p_AA_Aa_aa()[1]) * f2.p_AA_Aa_aa()[2]
                adder_6 = (Decimal(n_c.get(CLDCalculation.F_Aabb, 0)) / self.n - p_Aa_x_p_bb) ** 2 / p_Aa_x_p_bb

                p_aa_x_p_BB = Decimal(f1.p_AA_Aa_aa()[2]) * f2.p_AA_Aa_aa()[0]
                adder_7 = (Decimal(n_c.get(CLDCalculation.F_aaBB, 0)) / self.n - p_aa_x_p_BB) ** 2 / p_aa_x_p_BB

                p_aa_x_p_Bb = Decimal(f1.p_AA_Aa_aa()[2]) * f2.p_AA_Aa_aa()[1]
                adder_8 = (Decimal(n_c.get(CLDCalculation.F_aaBb, 0)) / self.n - p_aa_x_p_Bb) ** 2 / p_aa_x_p_Bb

                p_aa_x_p_bb = Decimal(f1.p_AA_Aa_aa()[2]) * f2.p_AA_Aa_aa()[2]
                adder_9 = (Decimal(n_c.get(CLDCalculation.F_aabb, 0)) / self.n - p_aa_x_p_bb) ** 2 / p_aa_x_p_bb

                cld2 = self.n * (
                    adder_1 + adder_2 + adder_3 + adder_4 + adder_5 + adder_6 + adder_7 + adder_8 + adder_9)

                if cld2 > cld2_max:
                    cld2_max = cld2
        return cld2_max

    def __n_AB(self, snp_feature1, snp_feature2):
        n_c = self.__n_count(snp_feature1, snp_feature2)
        return 2 * n_c.get(self.F_AABB, 0) + n_c.get(self.F_AABb, 0) + n_c.get(self.F_AaBB, 0) + Decimal(0.5) * n_c.get(
            self.F_AaBb, 0)

    def __n_Ab(self, snp_feature1, snp_feature2):
        n_c = self.__n_count(snp_feature1, snp_feature2)
        return 2 * n_c.get(self.F_AAbb, 0) + n_c.get(self.F_AABb, 0) + n_c.get(self.F_Aabb, 0) + Decimal(0.5) * n_c.get(
            self.F_AaBb, 0)

    def __n_aB(self, snp_feature1, snp_feature2):
        n_c = self.__n_count(snp_feature1, snp_feature2)
        return 2 * n_c.get(self.F_aaBB, 0) + n_c.get(self.F_aaBb, 0) + n_c.get(self.F_AaBB, 0) + Decimal(0.5) * n_c.get(
            self.F_AaBb, 0)

    def __n_ab(self, snp_feature1, snp_feature2):
        n_c = self.__n_count(snp_feature1, snp_feature2)
        return 2 * n_c.get(self.F_aabb, 0) + n_c.get(self.F_aaBb, 0) + n_c.get(self.F_Aabb, 0) + Decimal(0.5) * n_c.get(
            self.F_AaBb, 0)

    def __n_count(self, snp_feature1, snp_feature2):
        n_count = {}
        for i in range(self.n):
            key = (snp_feature1[i], snp_feature2[i])
            n_count[key] = n_count.get(key, 0) + 1
        return n_count
