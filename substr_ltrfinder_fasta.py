#!/usr/bin/python3

import os
import sys
import csv
import itertools
from time import time
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import Alphabet

if len(sys.argv) < 3:
    print('Usage: ', str(sys.argv[0]), '[LTR FINDER FILE] [FASTA FILE] [OUTPUT NAME]')
else:
    start_time = time()

    ltrfinder_file = str(sys.argv[1])
    fasta_file = str(sys.argv[2])

    records = SeqIO.index(fasta_file, 'fasta')
    output = open(sys.argv[3], 'w')

    print('Using', ltrfinder_file, 'LTR Finder and', fasta_file, 'FASTA File')

    with open(ltrfinder_file, 'r') as file:
        for line1, line2, line3 in itertools.zip_longest(*[file]*3):
            id = line1.split()[1]
            begin = line2.split()[2]
            end = line2.split()[4]

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

    print('\nCheck the results in ', sys.argv[3], '\n', sep='')

    output.close()
    records.close()

    end_time = time()

    print('Took %.3f seconds...\n' % (end_time - start_time))
