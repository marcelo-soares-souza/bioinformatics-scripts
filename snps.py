#!/usr/bin/env python3

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

    print('Using Prefix', prefix, '\n')

    for suffix, filename in filename.items():

        with open(filename, 'r') as csv_data:
            reader = csv.reader(csv_data, delimiter='\t')

            for data in reader:
                sequence = data[0]
                position = data[1]

                sequences[sequence][position].append(suffix)


    tools = defaultdict(list)

    for sequence, list in sequences.items():
        for key, suffix in list.items():
            tools[', '.join(suffix)].append(sequence + ' ' + key)

    tools_sorted = sorted(tools, key=len, reverse=True)

    output = open(prefix + '.result', 'w')

    for tool in tools_sorted:
        output.write('Found on %s\n\n' % (tool))

        for seq_pos in tools[tool]:
            output.write("%s\n" % (seq_pos))

        output.write('\n')

    output.close()

    end_time = time()


    print('\nTook %.3f seconds...\n' % (end_time - start_time))
