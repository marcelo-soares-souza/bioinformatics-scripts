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

if len(sys.argv) < 5:
    print('Usage:', str(sys.argv[0]),
          '[CSV FILE] [FASTQ FILE] [PIDENT] [QCOVS]')
else:
    start_time = time()

    csv_filename = str(sys.argv[1])
    fastq_filename = str(sys.argv[2])
    arg_pident = float(sys.argv[3])
    arg_qcovs = float(sys.argv[4])

    output_result_filename = os.path.splitext(sys.argv[2])[
        0] + '.clean.' + datetime.datetime.now().strftime('%Y%m%d') + '.fastq'
    output_removed_filename = os.path.splitext(sys.argv[2])[
        0] + '.removed.' + datetime.datetime.now().strftime('%Y%m%d') + '.fastq'

    output_removed = open(output_removed_filename, 'w')

    write = False
    removed_sequences = 0

    print('\nReading', fastq_filename, 'FASTQ File')
    records = SeqIO.to_dict(SeqIO.parse(fastq_filename, 'fastq'))

    print('\nUsing', csv_filename, 'CSV Table and', fastq_filename, 'FASTQ File')

    with open(csv_filename, 'r') as csv_file:
        reader = csv.reader(csv_file, delimiter='\t')

        for data in reader:
            id = data[0]
            subject = data[2]
            pident = float(data[11])
            qcovs = float(data[14])

            if pident >= arg_pident and qcovs >= arg_qcovs:
                if id in records.keys():

                    print('\nRemoving', id, 'with PIDENT',
                          pident, 'and QCOVS', qcovs)

                    header = '%s %s %s %s' % (
                        id, subject, str(pident), str(qcovs))
                    record = SeqRecord(
                        Seq(str(records[id].seq), Alphabet()), id=str(header), description='')
                    record.letter_annotations["phred_quality"] = records[
                        id].letter_annotations["phred_quality"]

                    SeqIO.write(record, output_removed, 'fastq')

                    del records[id]

                    write = True
                    removed_sequences = removed_sequences + 1

    output_removed.close()

    print('Writing Results...')

    if write:
        output = open(output_result_filename, 'w')
        [SeqIO.write(record, output, 'fastq')
         for (id, record) in records.items()]
        output.close()

        print('\nCheck the results in ', output_result_filename, '\n', sep='')
        print('Removed Sequences (', removed_sequences, ') in ',
              output_removed_filename, '\n', sep='')
    else:
        print('No sequences found\n')

    end_time = time()

    print('Took %.3f seconds...\n' % (end_time - start_time))
