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
        sys.argv[0]), '[CSV FILE] [FASTA FILE] [UPSTREAM/DOWNSTREAM BP VALUE]')
else:
    start_time = time()

    csv_file = str(sys.argv[1])
    fasta_file = str(sys.argv[2])
    up_down_value = int(sys.argv[3])

    records = SeqIO.index(fasta_file, 'fasta')
    output = open(os.path.splitext(sys.argv[2])[
                  0] + '.result.' + datetime.datetime.now().strftime('%Y%m%d') + '.fasta', 'w')

    print('Using', csv_file, 'CSV Table and', fasta_file,
          'FASTA File', 'With Up/Down Stream of', up_down_value)

    with open(csv_file, 'r') as csv_file:
        reader = csv.reader(csv_file)

        for name, id, begin, end in reader:
            end = int(begin)
            begin = int(begin) - up_down_value

            print('\nProcessing', name, id, 'of size', len(
                records[id]), 'starting in', begin, 'ending', end)

            if int(begin) < int(end):
                print('Normal Slicing...')
                seq = records[id].seq[int(begin):int(end)]
                id_name = id
            else:
                print('Reverse Complement...')
                seq = records[id].seq[int(end):int(begin)].reverse_complement()
                id_name = id + '_RC'

            header = '%s|size%s[%s_to_%s](%s nts)' % (
                id_name, len(records[id]), begin, end, len(seq))

            print(header)

            record = SeqRecord(Seq(str(seq), Alphabet()),
                               id=str(header), description='')

            SeqIO.write(record, output, 'fasta')

    print('\nCheck the results in ', os.path.splitext(
        sys.argv[2])[0], '.result.fasta\n', sep='')

    output.close()
    records.close()

    end_time = time()

    print('Took %.3f seconds...\n' % (end_time - start_time))
