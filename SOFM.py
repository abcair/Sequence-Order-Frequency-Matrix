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

If you use this library, please cite the following papaers:

1. Junjie Chen, Mingyue Guo, Xiaolong Wang, BinLiu*. Protein remote homology detection and
   fold recognition based on Sequence-Order Frequency Matrix[J]. IEEE/ACM Transactions on
   Computational Biology and Bioinformatics, 2017.
2. Junjie Chen, Mingyue Guo, Xiaolong Wang, BinLiu*. SOFM-Top: Protein Remote Homology
   Detection and Fold Recognition Based on Sequence-Order Frequency Matrix[C]//International
   Conference on Intelligent Computing. Springer, Cham, 2017: 469-480.

Last Motified Time: 09/08/2017 

"""

import os
import sys
import argparse
import subprocess
import ConfigParser
import itertools
import collections
import numpy as np
import xml.etree.ElementTree as ET
from Bio import SeqIO

# 20 standard amino acids
alphabet = 'ACDEFGHIKLMNPQRSTVWY'

def Kmer2Index(kmer):
    kmer_tuple = []
    kmer_tuple = list(itertools.product(alphabet, repeat=kmer))
    kmer_string = [''.join(x) for x in kmer_tuple]
    kmer2index = {kmer_string[i]: i for i in range(0, len(kmer_string))}
    return kmer2index


def run_search(FASTAfile):
    # temp file names
    xml_file = os.path.join(filepath, shortname + '.xml')

    # run psi-blast command
    print "Running PSIBLAST search..."
    outfmt_type = 5
    num_iter = 3
    evalue_threshold = 0.001
    threads = 5    
    cmd = ' '.join([PSIBLAST,
                    '-query ' + FASTAfile,
                    '-db ' + BLAST_DB,
                    '-out ' + xml_file,
                    '-evalue ' + str(evalue_threshold),
                    '-num_iterations ' + str(num_iter),
                    '-outfmt ' + str(outfmt_type),
                    '-num_threads ' + str(threads)]
                   )    
    return_code = subprocess.call(cmd, shell=True)
    
    print 'Parsing xml file...'
    # parser the xml file output by PSIBLAST
    tree = ET.ElementTree(file=xml_file)
    query_def = tree.find('BlastOutput_query-def').text
    query_len = tree.find('BlastOutput_query-len').text
    # get the last iteration
    iteration = tree.findall('BlastOutput_iterations/Iteration')[-1]
    Iteration_hits = iteration.find('Iteration_hits')

    # store searched homology proteins
    MSA = []
    for Hit in Iteration_hits.getchildren():
        Hsp_evalue = Hit.find('Hit_hsps/Hsp/Hsp_evalue').text

        # only parser the hits that e-value < threshold
        if float(Hsp_evalue) > evalue_threshold:
            continue

        Hit_num = Hit.find('Hit_num').text
        Hit_id = Hit.find('Hit_id').text
        Hit_def = Hit.find('Hit_def').text

        Hsp_query_from = Hit.find('Hit_hsps/Hsp/Hsp_query-from').text
        Hsp_query_to = Hit.find('Hit_hsps/Hsp/Hsp_query-to').text
        Hsp_hit_from = Hit.find('Hit_hsps/Hsp/Hsp_hit-from').text
        Hsp_hit_to = Hit.find('Hit_hsps/Hsp/Hsp_hit-to').text
        Hsp_qseq = Hit.find('Hit_hsps/Hsp/Hsp_qseq').text
        Hsp_hseq = Hit.find('Hit_hsps/Hsp/Hsp_hseq').text

        # alignment sequence by add prefix, suffix
        prefix = "-" * (int(Hsp_query_from) - 1)
        suffix = "-" * (int(query_len) - int(Hsp_query_to))

        # delete the space in protein_name and the corresponding position of
        # hits
        pos = -1
        for aa in Hsp_qseq:
            pos = pos + 1
            if aa == '-':
                Hsp_hseq = Hsp_hseq[:pos] + '*' + Hsp_hseq[pos + 1:]
        Hsp_hseq = Hsp_hseq.replace('*', '')

        if 'X' in Hsp_hseq:
            Hsp_hseq = Hsp_hseq.replace('X', '-')
        if 'B' in Hsp_hseq:
            Hsp_hseq = Hsp_hseq.replace('B', '-')
        if 'Z' in Hsp_hseq:
            Hsp_hseq = Hsp_hseq.replace('Z', '-')
        if 'U' in Hsp_hseq:
            Hsp_hseq = Hsp_hseq.replace('U', '-')
        if 'J' in Hsp_hseq:
            Hsp_hseq = Hsp_hseq.replace('J', '-')
        if 'O' in Hsp_hseq:
            Hsp_hseq = Hsp_hseq.replace('O', '-')

        # combine prefix, modified hits, suffix
        hit_sequence = prefix + Hsp_hseq + suffix
        # print hit_sequence

        # append in MSA
        MSA.append(hit_sequence)

    if MSA == []:
        # append the protein-self
        ff = open(fasta_file, 'r')
        ff.readline()  # skip the id
        fasta_seq = ff.readline().strip().upper()
        ff.close()

        if 'X' in fasta_seq:
            fasta_seq = fasta_seq.replace('X', '-')
        if 'B' in fasta_seq:
            fasta_seq = fasta_seq.replace('B', '-')
        if 'Z' in fasta_seq:
            fasta_seq = fasta_seq.replace('Z', '-')
        if 'U' in fasta_seq:
            fasta_seq = fasta_seq.replace('U', '-')
        if 'J' in fasta_seq:
            fasta_seq = fasta_seq.replace('J', '-')
        if 'O' in fasta_seq:
            fasta_seq = fasta_seq.replace('O', '-')

        MSA.append(fasta_seq)

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


def WriteSOFM(FASTAseq, SOFMfile, SOFM, kmer2index):
    headers = sorted(kmer2index.iteritems(), key=lambda d: d[1])

    with open(SOFMfile, 'w') as f:
        f.write("{:<9s}".format(''))
        for item in headers:
            f.write("{:<8s}".format(item[0]))
        f.write('\n')

        for i in xrange(0, SOFM.shape[0]):
            # rower number and fasta sequence
            f.write("{:<4s}".format(str(i + 1)) + "{:<5s}".format(FASTAseq[i]))
            for j in xrange(0, SOFM.shape[1]):
                f.write("{:<8s}".format(str(SOFM[i, j])))
            f.write('\n')


def main(FASTAfile, SOFMfile, kmer):
    # a dictionary that contains all amino acid substrings with length k
    kmer2index = Kmer2Index(kmer)

    
    MSA = run_search(FASTAfile)
    
    print 'Generating SOFM...'
    if MSA != []:
        SOFM = GenerateSOFM(MSA, kmer, kmer2index)

    print 'Writing output file...'
    # extract fasta sequence
    FASTAseq = list(SeqIO.parse(open(FASTAfile), 'fasta'))[0].seq.upper()
    WriteSOFM(FASTAseq, SOFMfile, SOFM, kmer2index)


if __name__ == '__main__':
    # parameters
    parser = argparse.ArgumentParser(
        description="Sequence-Order Frequency Matrix (SOFM) is a novel protein profile, which can achieve more information content than traditional profiles.")
    parser.add_argument('-i', '--input', required=True,
                        help='input a single protein sequence in FASTA format')
    parser.add_argument('-o', '--output', default=None, required=False,
                        help='output SOFM file, default filename: {input}.sofm{k}')
    parser.add_argument('-k', '--kmer', type=int, default=3, required=False,
                        help='length of substrings, default value: k=3')
    parser.add_argument('--version', action='version', version='%(prog)s 1.0')
    args = parser.parse_args()

    FASTAfile = args.input
    SOFMfile = args.output
    kmer = args.kmer

    (filepath, filename) = os.path.split(FASTAfile)
    (shortname, extension) = os.path.splitext(filename)
    if not SOFMfile:
        SOFMfile = os.path.join(filepath, shortname + '.sofm' + str(kmer))

    conf = ConfigParser.ConfigParser()
    conf.read('PATH.conf')
    PSIBLAST = os.path.join(conf.get('PSIBLAST', 'NCBI_BLAST_BIN'),
                            'psiblast')
    BLAST_DB = conf.get('PSIBLAST', 'BLAST_DATABASE')

    main(FASTAfile, SOFMfile, kmer)
    print 'DONE!'
