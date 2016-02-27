from __future__ import print_function
from collections import OrderedDict

import os

import numpy
import re

from snp import SnpFeature

map_matrix = numpy.genfromtxt('test.map', str)
ped_matrix = numpy.genfromtxt('test.ped', str)

map_matrix = map_matrix[:, -1]
ped_matrix = ped_matrix[:, 6:]

GROUP_RANGE = 100000
SNP_SIZE = 2


def group_and_save_snps():
    __remove_output()
    for i in range(0, int(len(ped_matrix[0]) / SNP_SIZE)):
        snp_i = numpy.transpose(ped_matrix[:, i * SNP_SIZE:(i + 1) * SNP_SIZE])
        group_num = int(map_matrix[i])
        group_id = int(int(group_num) / GROUP_RANGE)
        __save_snp(group_id, snp_i)

        # fixme for debug only >>>

        snp_feature = SnpFeature(snp_i)
        print(snp_feature.p_A())
        print(snp_feature.p_AA_Aa_aa())
        # fixme for debug only <<<

    return __get_output_file_names()


def __remove_output():
    for f in os.listdir('.'):
        if re.search('.+\.out$', f):
            os.remove(f)


def __save_snp(snp_id, snp):
    with open('group.' + str(snp_id) + '.out', 'a') as f:
        f.write(str(snp))
        f.write('\n\n')


def __get_output_file_names():
    output_map = {}
    for file_name in os.listdir('.'):
        if re.search('.+\.out$', file_name):
            file_id = int(file_name[len('group.'): -len('.out')])
            output_map[file_id] = file_name
    return output_map


def main_process():
    output_map = group_and_save_snps()
    output_map.popitem()


if __name__ == '__main__':
    main_process()
