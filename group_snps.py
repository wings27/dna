from __future__ import print_function

import numpy

map_matrix = numpy.genfromtxt('test.map', str)
ped_matrix = numpy.genfromtxt('test.ped', str)

map_matrix = map_matrix[:, -1]
ped_matrix = ped_matrix[:, 6:]

GROUP_RANGE = 100000
SNP_SIZE = 2


def save_snp(id, snp):
    with open('group.' + str(id) + '.out', 'a') as f:
        f.write(str(snp))
        f.write('\n\n')


def pA(snp):
    snp00 = snp[0][0]
    if snp00 == 'A' or snp00 == 'G':
        count_num = count_item_in_snp(snp, 'A')
    elif snp00 == 'T':
        count_num = count_item_in_snp(snp, 'C')
    else:
        count_a = count_item_in_snp(snp, 'A')
        count_num = count_a if count_a != 0 else count_item_in_snp(snp, 'C')
    return float(count_num) / (len(snp[0]) * SNP_SIZE)


def pAA_Aa_aa(snp):
    snp00 = snp[0][0]
    if snp00 == 'A' or snp00 == 'G':
        result = count_pair_in_snp(snp, 'A', 'G')
    elif snp00 == 'T':
        result = count_pair_in_snp(snp, 'C', 'T')
    else:
        count_a = count_item_in_snp(snp, 'A')
        result = count_pair_in_snp(snp, 'A', 'C') if count_a != 0 else count_pair_in_snp(snp, 'C', 'T')

    return tuple(map(lambda x: float(x) / len(snp[0]), result))


def count_item_in_snp(snp, item):
    return list(snp[0]).count(item) + list(snp[1]).count(item)


def count_pair_in_snp(snp, large_item, small_item):
    l_key = large_item * 2
    ls_key = large_item + small_item
    sl_key = small_item + large_item
    s_key = small_item * 2
    count_dict = {l_key: 0, ls_key: 0, sl_key: 0, s_key: 0}
    for i in range(0, len(snp[0])):
        symbol = snp[0][i] + snp[1][i]
        count_dict[symbol] += 1

    return count_dict[l_key], count_dict[ls_key] + count_dict[sl_key], count_dict[s_key]


def extract_feature(snp):
    snp00 = snp[0][0]
    if snp00 == 'A' or snp00 == 'G':
        result = extract_feature_ls(snp, 'A', 'G')
    elif snp00 == 'T':
        result = extract_feature_ls(snp, 'C', 'T')
    else:
        count_a = count_item_in_snp(snp, 'A')
        result = extract_feature_ls(snp, 'A', 'C') if count_a != 0 else extract_feature_ls(snp, 'C', 'T')
    return result


def extract_feature_ls(snp, large_item, small_item):
    return [extract_feature_from_symbol(snp[0][i] + snp[1][i], large_item, small_item) for i in range(len(snp[0]))]


def extract_feature_from_symbol(symbol, large_item, small_item):
    if symbol == large_item * 2:
        return 1
    elif symbol == small_item * 2:
        return -1
    else:
        return 0


def main_process():
    for i in range(0, int(len(ped_matrix[0]) / SNP_SIZE)):
        snp_i = numpy.transpose(ped_matrix[:, i * SNP_SIZE:(i + 1) * SNP_SIZE])
        group_num = int(map_matrix[i])
        group_id = int(int(group_num) / GROUP_RANGE)
        save_snp(group_id, snp_i)
        print(extract_feature(snp_i))


if __name__ == '__main__':
    main_process()
