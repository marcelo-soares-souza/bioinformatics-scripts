#!/usr/bin/python3

import os
import sys
import csv
from time import time
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import Alphabet

if len(sys.argv) < 2:
    print('Usage: ', str(sys.argv[0]), '[CSV FILE] [FASTA FILE]')
else:
    start_time = time()

    csv_file = str(sys.argv[1])
    fasta_file = str(sys.argv[2])

    records = SeqIO.index(fasta_file, 'fasta')
    output = open(os.path.splitext(sys.argv[2])[0] + '.result.fasta', 'w')

    print('Using', csv_file, 'CSV Table and', fasta_file, 'FASTA File')

    with open(csv_file, 'r') as csv_file:
        reader = csv.reader(csv_file)

        for id, begin, end in reader:
            print('\nProcessing', id, 'of size', len(records[id]), 'starting in', begin, 'ending', end)

            if int(begin) < int(end):
                print('Normal Slicing...')
                seq = records[id].seq[int(begin):int(end)]
            else:
                print('Reverse Complement...')
                seq = records[id].seq[int(end):int(begin)].reverse_complement()

            header = '%s|size%s[%s_to_%s](%s nts)' % (id, len(records[id]), begin, end, len(seq))

            record = SeqRecord(Seq(str(seq), Alphabet()), id=str(header), description='')

            SeqIO.write(record, output, 'fasta')

    print('\nCheck the results in ', os.path.splitext(sys.argv[2])[0], '.result.fasta\n', sep='')

    output.close()
    records.close()

    end_time = time()

    print('Took %.3f seconds...\n' % (end_time - start_time))
