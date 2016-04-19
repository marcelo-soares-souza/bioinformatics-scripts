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

if len(sys.argv) < 3:
    print('Usage:', str(sys.argv[0]), '[CSV FILE] [FASTA FILE]')
else:
    start_time = time()

    csv_filename = str(sys.argv[1])
    fasta_filename = str(sys.argv[2])
    output_result_filename = os.path.splitext(sys.argv[2])[0] + '.result.' + datetime.datetime.now().strftime('%Y%m%d') + '.fasta'
    output_removed_filename = os.path.splitext(sys.argv[2])[0] + '.removed.' + datetime.datetime.now().strftime('%Y%m%d') + '.fasta'

    output_removed = open(output_removed_filename, 'w')

    write = False
    removed_sequences = 0

    print('\nReading', fasta_filename, 'FASTA File')
    records = SeqIO.to_dict(SeqIO.parse(fasta_filename, 'fasta'))

    print('\nUsing', csv_filename, 'CSV Table and', fasta_filename, 'FASTA File')

    with open(csv_filename, 'r') as csv_file:
        reader = csv.reader(csv_file)

        for id in reader:
            id = id[0]
            if id in records.keys():
              print('\nRemoving', id)

              SeqIO.write(records[id], output_removed, 'fasta')

              del records[id]

              write = True
              removed_sequences = removed_sequences + 1

    if write:
        output = open(output_result_filename, 'w')
        [ SeqIO.write(record, output, 'fasta') for (id, record) in records.items() ]
        output.close()
        print('\nCheck the results in ', output_result_filename, '\n', sep='')
        print('\nRemoved Sequences (',removed_sequences,') in ', output_removed_filename, '\n', sep='')
    else:
        print('No sequences found\n')
    
    end_time = time()

    print('Took %.3f seconds...\n' % (end_time - start_time))

    output_removed.close()