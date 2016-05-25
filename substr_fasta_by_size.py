#!/usr/bin/python3

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
    print('Usage:', str(sys.argv[0]), '[SIZE] [FASTA FILE] [OUTPUT]')
else:
    start_time = time()

    size = int(sys.argv[1])
    fasta_file = str(sys.argv[2])

    records = SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta'))
    output = open(sys.argv[3], "w")

    print('Using size', size, 'and', fasta_file, 'FASTA File')

    for id in records.keys():
        print('\nRecord', records[id])

        seq = records[id].seq[:size]
        header = '%s' % (id)

        record = SeqRecord(Seq(str(seq), Alphabet()), id=str(header), description='')

        SeqIO.write(record, output, 'fasta')

    print('\nCheck the results in ', sys.argv[3], '\n', sep='')

    output.close()

    end_time = time()

    print('Took %.3f seconds...\n' % (end_time - start_time))
