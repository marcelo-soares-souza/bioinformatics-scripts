#!/usr/bin/python3

# (C) 2016 Marcelo Soares Souza <marcelo.soares@colaborador.embrapa.br>
# This program is licensed under a LGPLv3 License.

import os
import sys
import csv
import datetime
from collections import defaultdict
from time import time

def getKey(item):
    return item[1]

if len(sys.argv) < 2:
    print('Usage:', str(sys.argv[0]), '[PREFIX FILENAME]')
else:
    start_time = time()

    prefix = str(sys.argv[1])

    filename = {}
    filename['bowtie_freebayes'] = prefix + '_' + 'bowtie_freebayes.csv'
    filename['bowtie_GATK'] = prefix + '_' + 'bowtie_GATK.csv'
    filename['bowtie_samtools'] = prefix + '_' + 'bowtie_samtools.csv'
    filename['bwa_freebayes'] = prefix + '_' + 'bwa_freebayes.csv'
    filename['bwa_GATK'] = prefix + '_' + 'bwa_GATK.csv'
    filename['bwa_samtools'] = prefix + '_' + 'bwa_samtools.csv'

    sequences = defaultdict(lambda: defaultdict(list))
    # sequences = defaultdict(list)

    print('Using Prefix', prefix, '\n')

    for suffix, filename in filename.items():

        with open(filename, 'r') as csv_data:
            reader = csv.reader(csv_data, delimiter='\t')

            for data in reader:
                sequence = data[0]
                position = data[1]

                sequences[sequence][position].append(suffix)


    s = defaultdict(list)

    for sequence, list in sequences.items():
        for key, suffix in list.items():
            s[', '.join(suffix)].append(sequence + ' ' + key)

    ss = sorted(s, key=len, reverse=True)

    for x in ss:
        print('\n')
        print(x)
        print('\n')
        for a in s[x]:
          print(a)
        print('\n')

    end_time = time()

    print('\nTook %.3f seconds...\n' % (end_time - start_time))
