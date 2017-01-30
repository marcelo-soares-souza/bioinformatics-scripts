#!/usr/bin/env python3

# (C) 2016 Marcelo Soares Souza <marcelo@riseup.net>
# This program is licensed under a LGPLv3 License.

import os
import sys
import csv
import datetime
from time import time
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import Alphabet

if len(sys.argv) < 4:
    print('Usage:', str(
        sys.argv[0]), '[CSV FILE] [FASTA FILE] [DOWNSTREAM BP VALUE]')
else:
    start_time = time()

    csv_file = str(sys.argv[1])
    fasta_file = str(sys.argv[2])
    down_value = int(sys.argv[3])

    records = SeqIO.index(fasta_file, 'fasta')
    output = open(os.path.splitext(sys.argv[2])[0] + '.result.' + datetime.datetime.now().strftime('%Y%m%d') + '.fasta', 'w')

    print('Using', csv_file, 'CSV Table and', fasta_file, 'FASTA File', 'With Up/Down Stream of', down_value)

    with open(csv_file, 'r') as csv_file:

        reader = csv.reader(csv_file, delimiter=';')

        for name, id, begin in reader:
            end = int(begin)
            begin = int(begin) - int(down_value)

            print('\nProcessing', name, id, 'of size', len(records[id]), 'starting in', begin, 'ending', end)

            print('Normal Slicing...')
            seq = records[id].seq[int(begin):int(end)]
            id_name = name + '_' + id

            # header = '%s|size%s[%s_to_%s](%s nts)' % (id_name, len(records[id]), begin, end, len(seq))
            header = '%s' % (id_name)

            print(header)

            record = SeqRecord(Seq(str(seq), Alphabet()), id=str(header), description='')

            SeqIO.write(record, output, 'fasta')

    print('\nCheck the results in ', os.path.splitext(sys.argv[2])[0], '.result.fasta\n', sep='')

    output.close()
    records.close()

    end_time = time()

    print('Took %.3f seconds...\n' % (end_time - start_time))
