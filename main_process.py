from __future__ import print_function

import os
import re

import itertools

from snp import *


class MatrixContainer:
    map_matrix = numpy.genfromtxt('newtest.map', str)
    ped_matrix = numpy.genfromtxt('newtest.ped', str)

    map_matrix = map_matrix[:, -1]
    ped_matrix = ped_matrix[:, 6:]


GROUP_RANGE = 100000
SNP_SIZE = 2


def group_and_save_feature():
    __remove_output()
    prefix = 'feature.'

    p_matrix = MatrixContainer.ped_matrix
    m_matrix = MatrixContainer.map_matrix
    for i in range(0, int(len(p_matrix))):
        if i >= len(p_matrix):
            break
        l = list(p_matrix[i])
        if '0' in l:
            p_matrix = numpy.delete(p_matrix, i, 0)

    for i in range(0, int(len(p_matrix[0]) / SNP_SIZE)):
        snp_i = numpy.transpose(p_matrix[:, i * SNP_SIZE:(i + 1) * SNP_SIZE])
        group_num = int(m_matrix[i])
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

    smaller_size = min(len(output_map_keys), 20)
    result_matrix = numpy.zeros((smaller_size, smaller_size))

    for cbn in itertools.combinations(range(smaller_size), 2):
        x = cbn[0]
        y = cbn[1]

        file1 = output_map[output_map_keys[x]]
        file2 = output_map[output_map_keys[y]]

        features1 = FileHelper.load_feature_group(file1)
        features2 = FileHelper.load_feature_group(file2)
        cld_cal = CLDCalculation(features1, features2)
        cld_1 = cld_cal.cld_1()

        result_matrix[x][y] = cld_1
        result_matrix[y][x] = cld_1

    for diagonal in range(smaller_size):
        file = output_map[output_map_keys[diagonal]]
        features = FileHelper.load_feature_group(file)
        cld_cal = CLDCalculation(features, features)
        cld_1 = cld_cal.cld_1()
        result_matrix[diagonal][diagonal] = cld_1

    print(result_matrix)


if __name__ == '__main__':
    main_process()
