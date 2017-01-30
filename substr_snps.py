#!/usr/bin/env python3

# (C) 2016 Marcelo Soares Souza <marcelo@riseup.net>
# This program is licensed under a LGPLv3 License.

import os
import sys
import csv
import datetime
from time import time
from Bio import SeqIO

if len(sys.argv) < 4:
    print('Usage:', str(
        sys.argv[0]), '[ORGANISM] [CSV FILE] [FASTA FILE] [UPSTREAM/DOWNSTREAM BP VALUE]')
else:
    start_time = time()

    organism = str(sys.argv[1])
    csv_file = str(sys.argv[2])
    fasta_file = str(sys.argv[3])
    up_down_value = int(sys.argv[4])

    records = SeqIO.index(fasta_file, 'fasta')
    output = open(os.path.splitext(sys.argv[2])[0] + '.result.' + datetime.datetime.now().strftime('%Y%m%d') + '.snps', 'w')

    print('Using', csv_file, 'CSV Table and', fasta_file, 'FASTA File', 'With Up/Down Stream of', up_down_value)

    with open(csv_file, 'r') as csv_file:
        reader = csv.reader(csv_file, delimiter='\t')

        for chr, pos, ref, alt in reader:
            begin = int(pos) - up_down_value
            end = int(pos) + up_down_value

            seq = str(records[chr].seq[int(begin):int(end)])
            seq = seq[:up_down_value] + '[' + ref + '/' + alt + ']' + seq[up_down_value:]
            snpid = '%s_%s_%s' % (organism, chr, pos)
            info = '%s;%s;%s;%s\n' % (snpid, seq, chr, pos)

            output.write(info)

    print('\nCheck the results in ', os.path.splitext(sys.argv[2])[0], '.result.csv\n', sep='')

    output.close()
    records.close()

    end_time = time()

    print('Took %.3f seconds...\n' % (end_time - start_time))
