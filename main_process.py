from __future__ import print_function

import os
import re

import numpy

from snp import *

map_matrix = numpy.genfromtxt('test.map', str)
ped_matrix = numpy.genfromtxt('test.ped', str)

map_matrix = map_matrix[:, -1]
ped_matrix = ped_matrix[:, 6:]

GROUP_RANGE = 100000
SNP_SIZE = 2


def group_and_save_feature():
    __remove_output()
    prefix = 'feature.'
    for i in range(0, int(len(ped_matrix[0]) / SNP_SIZE)):
        snp_i = numpy.transpose(ped_matrix[:, i * SNP_SIZE:(i + 1) * SNP_SIZE])
        group_num = int(map_matrix[i])
        group_id = int(int(group_num) / GROUP_RANGE)

        snp_feature = SnpFeature.extract_feature(snp_i)

        __save_array(prefix + str(group_id), snp_feature)

    output_map = {}
    for file_name in os.listdir('.'):
        if re.search('.+\.out$', file_name):
            file_id = int(file_name[len(prefix): -len('.out')])
            output_map[file_id] = file_name
    return output_map


def __remove_output():
    for f in os.listdir('.'):
        if re.search('.+\.out$', f):
            os.remove(f)


def __save_array(file_name, array):
    with open(file_name + '.out', 'a') as f:
        f.write(str(array)[1:-1])
        f.write('\n')


def main_process():
    output_map = group_and_save_feature()
    output_map_keys = sorted(output_map.keys())
    file1 = output_map[output_map_keys[0]]
    file2 = output_map[output_map_keys[1]]

    print(file1, file2)
    features1 = FileHelper.load_feature_group(file1)
    features2 = FileHelper.load_feature_group(file2)
    cld_cal = CLDCalculation(features1, features2)
    print(cld_cal.cld_1())


if __name__ == '__main__':
    main_process()
