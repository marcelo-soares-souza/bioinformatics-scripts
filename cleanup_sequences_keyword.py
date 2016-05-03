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
    print('Usage:', str(sys.argv[0]), '[FASTQ FILE] [KEYWORD]')
else:
    start_time = time()

    fastq_filename = str(sys.argv[1])
    keyword = str(sys.argv[2])

    output_result_filename = os.path.splitext(sys.argv[1])[0] + '.without_' + keyword + '.fastq'

    print('\nReading', fastq_filename, 'FASTQ File')

    records = SeqIO.to_dict(SeqIO.parse(fastq_filename, 'fastq'))

    for record in list(records.items()):
        if keyword.lower() in record[1].description.lower():
            print('\nRemoving', record[1].description)
            del records[record[0]]

    output = open(output_result_filename, 'w')
    [SeqIO.write(record, output, 'fastq') for (id, record) in records.items()]
    output.close()

    print('\nCheck the results in ', output_result_filename, '\n', sep='')

    end_time = time()

    print('\nTook %.3f seconds...\n' % (end_time - start_time))
