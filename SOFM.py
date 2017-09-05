"""
Copyright (C) 2016-2017 Junjie Chen (junjie.chen.hit@gmail.com),

Harbin Institute of Technology, Shenzhen, China
Bioinformatics Group.
The software is maintained by Junjie Chen.

This program is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation; either version 2 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for more details.
 
"""


import sys
import argparse

import itertools
import collections
import numpy as np


def Kmer2Index(alphabet, kmer):
    kmer_tuple = []
    kmer_tuple = list(itertools.product(alphabet, repeat=kmer))
    kmer_string = [''.join(x) for x in kmer_tuple]
    kmer2index = {kmer_string[i]: i for i in range(0, len(kmer_string))}
    return kmer2index


def readMSA(msa_file):
    MSA = []
    with open(msa_file) as f:
        for line in f:
            MSA.append(line.strip())

    return MSA


def GenerateSOFM(MSA, kmer, kmer2index):
    # initialize SOFM    
    row_size = len(kmer2index)
    column_size = len(MSA[0]) - kmer + 1
    SOFM = np.zeros((row_size, column_size))

    # FREQUENCY MATRIX
    for col in xrange(column_size):
        substrings_at_a_position = []
        for row in MSA:
            substring = row[col:col + kmer]
            if '-' not in substring:
                substrings_at_a_position.append(substring)

        # count frequency
        substring2frequency = collections.Counter(substrings_at_a_position)

        for substring in substring2frequency:
            SOFM[kmer2index[substring], col] = substring2frequency[substring]

    return SOFM.transpose()


def WriteSOFM(SOFMfile, SOFM, kmer2index):
    headers = sorted(kmer2index.iteritems(), key=lambda d: d[1])

    with open(SOFMfile, 'w') as f:
        f.write("{:<9s}".format(''))
        for item in headers:
            f.write("{:<8s}".format(item[0]))
        f.write('\n')

        for i in xrange(0, SOFM.shape[0]):
            f.write("{:<9s}".format(str(i + 1) + '   '))  # rower number
            for j in xrange(0, SOFM.shape[1]):
                f.write("{:<8s}".format(str(SOFM[i, j])))
            f.write('\n')


def main(MSAfile, SOFMfile, kmer, kmer2index):
    MSA = readMSA(MSAfile)
    if MSA != []:
        SOFM = GenerateSOFM(MSA, kmer, kmer2index)

    WriteSOFM(SOFMfile, SOFM, kmer2index)


if __name__ == '__main__':
    # parameters
    parser = argparse.ArgumentParser(description="program description")
    parser.add_argument('MSAfile', help='input MSA file')
    parser.add_argument('SOFMfile', help='output SOFM file')
    parser.add_argument('-k', '--kmer', type=int, default=3,
                        help='length of substrings')
    #parser.add_argument('-norm', type=bool, default=True, help='normlize the column')
    args = parser.parse_args()

    MSAfile = args.MSAfile
    SOFMfile = args.SOFMfile
    kmer = args.kmer
    # print MSAfile, SOFMfile, kmer

    alphabet = 'ACDEFGHIKLMNPQRSTVWY'
    # a dictionary that contains all amino acid substrings with length k
    kmer2index = Kmer2Index(alphabet, kmer)

    main(MSAfile, SOFMfile, kmer, kmer2index)
