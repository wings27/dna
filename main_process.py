from __future__ import print_function

import itertools
import os
import re

from matplotlib import pyplot

from snp import *


class MatrixContainer:
    @classmethod
    def ped_matrix(cls, directory):
        ped_matrix = numpy.genfromtxt(directory + '.ped', str)
        ped_matrix = ped_matrix[:, 6:]
        return ped_matrix

    @classmethod
    def map_matrix(cls, directory):
        map_matrix = numpy.genfromtxt(directory + '.map', str)
        map_matrix = map_matrix[:, -1]
        return map_matrix


CONFIG = {
    'DIRECTORY': 'newtest',
    'GROUP_RANGE': 100000,
    'SNP_SIZE': 2
}


def group_and_save_feature():
    directory = CONFIG['DIRECTORY']

    __remove_output(directory)
    __remove_output('.')

    if not os.path.exists(directory):
        os.makedirs(directory)

    prefix = 'feature.'

    m_matrix = MatrixContainer.map_matrix(directory)
    p_matrix = MatrixContainer.ped_matrix(directory)

    for i in range(0, int(len(p_matrix))):
        if i >= len(p_matrix):
            break
        l = list(p_matrix[i])
        if '0' in l:
            p_matrix = numpy.delete(p_matrix, i, 0)

    snp_size, group_range = CONFIG['SNP_SIZE'], CONFIG['GROUP_RANGE']
    for i in range(0, int(len(p_matrix[0]) / snp_size)):
        snp_i = numpy.transpose(p_matrix[:, i * snp_size:(i + 1) * snp_size])
        group_num = int(m_matrix[i])
        group_id = int(int(group_num) / group_range)

        snp_feature = SnpFeature.extract_feature(snp_i)

        feature_matrix = numpy.matrix(snp_feature)
        __save_array(directory + '/' + prefix + str(group_id), feature_matrix, '%d')

    output_map = {}
    for file_name in os.listdir(directory):
        if re.search('.+\.out$', file_name):
            file_id = int(file_name[len(prefix): -len('.out')])
            output_map[file_id] = directory + '/' + file_name
    return output_map


def __remove_output(directory):
    if not os.path.exists(directory):
        return
    for f in os.listdir(directory):
        if re.search('.+\.out$', f):
            os.remove(directory + '/' + f)


def __save_array(file_name, array, fmt):
    with open(file_name + '.out', 'ab') as f:
        numpy.savetxt(f, array, fmt, ',')


def __render_array(result_matrix, interpolation):
    ax = pyplot.figure().gca()
    img = ax.imshow(result_matrix, interpolation=interpolation)
    pyplot.colorbar(img)
    pyplot.show()


def main_process():
    output_map = group_and_save_feature()
    output_map_keys = sorted(output_map.keys())

    smaller_size = min(len(output_map_keys), 20000)
    result_matrix = numpy.zeros((smaller_size, smaller_size))

    for cbn in itertools.combinations(range(smaller_size), 2):
        x = cbn[0]
        y = cbn[1]

        file1 = output_map[output_map_keys[x]]
        file2 = output_map[output_map_keys[y]]

        features1 = FileHelper.load_feature_group(file1)
        features2 = FileHelper.load_feature_group(file2)
        cld_cal = CLDCalculation(features1, features2)
        # cld_1 = cld_cal.cld_1()

        #  result_matrix[x][y] = cld_1
        #  result_matrix[y][x] = cld_1

        # cld_2 = cld_cal.cld_2()

        # result_matrix[x][y] = cld_2
        # result_matrix[y][x] = cld_2

        temp = cld_cal.temp_AB()

        result_matrix[x][y] = temp
        result_matrix[y][x] = temp
    print(result_matrix)

    __save_array('n_AB', result_matrix, '%.8f')


if __name__ == '__main__':
    main_process()
