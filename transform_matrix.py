import os

import itertools
import numpy
import re

from snp import FileHelper

FILE_NAME = 'FILE_NAME'
FEATURE_OUTPUT_DIR = 'FEATURE_OUTPUT_DIR'
MATRIX_OUTPUT = 'MATRIX_OUTPUT'

CONFIG = {
    FILE_NAME: 'n_AB.out',
    FEATURE_OUTPUT_DIR: 'newtest',
    MATRIX_OUTPUT: 'new_n_AB',
}


def load_feature_list():
    prefix = 'feature.'
    feature_output_dir = CONFIG[FEATURE_OUTPUT_DIR]
    return [int(f[len(prefix): -len('.out')])
            for f in os.listdir(feature_output_dir)
            if re.search('.+\.out$', f)]


def transform_matrix():
    feature_list = load_feature_list()
    file_name = CONFIG[FILE_NAME]
    origin_matrix = numpy.loadtxt(file_name, delimiter=',')
    min_size = min(origin_matrix.shape[0], origin_matrix.shape[1])
    combination = itertools.combinations(range(min_size), 2)
    counter = 0

    columns = 3
    result_matrix = numpy.zeros((min_size * (min_size - 1) / 2, columns))
    for cbn in combination:
        x = cbn[0]
        y = cbn[1]
        key0 = feature_list[x]
        key1 = feature_list[y]
        value = origin_matrix[x][y]

        result_matrix[counter][0] = key0
        result_matrix[counter][1] = key1
        result_matrix[counter][2] = value
        counter += 1

    FileHelper.save_array(CONFIG[MATRIX_OUTPUT], result_matrix, '%.8f')


def matrix_intersection(file_name_0, file_name_1):
    matrix_0 = numpy.loadtxt(file_name_0, delimiter=',')
    matrix_1 = numpy.loadtxt(file_name_1, delimiter=',')
    # todo
    # matrix_0[]
    set_0 = set()



if __name__ == '__main__':
    transform_matrix()
    matrix_intersection(CONFIG[MATRIX_OUTPUT], 'chr22_22.100k.OEnorm')
