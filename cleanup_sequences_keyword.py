#!/usr/bin/python3

# (C) 2016 Marcelo Soares Souza <marcelo.soares@colaborador.embrapa.br>
# This program is licensed under a LGPLv3 License.

import os
import sys
from time import time
from Bio import SeqIO

if len(sys.argv) < 3:
    print('Usage:', str(sys.argv[0]), '[FASTQ FILE] [KEYWORD]')
else:
    start_time = time()

    fastq_filename = str(sys.argv[1])
    keyword = str(sys.argv[2])

    filename = {}
    filename['without_keyword'] = os.path.splitext(
        sys.argv[1])[0] + '.without_' + keyword + '.fastq'
    filename['with_keyword'] = os.path.splitext(
        sys.argv[1])[0] + '.with_' + keyword + '.fastq'

    print('\nReading', fastq_filename, 'FASTQ File')

    records = SeqIO.to_dict(SeqIO.parse(fastq_filename, 'fastq'))

    output = {}
    output['with_keyword'] = open(filename['with_keyword'], 'w')
    for record in list(records.items()):
        if keyword.lower() in record[1].description.lower():
            print('\nRemoving', record[1].description)

            SeqIO.write(record[1], output['with_keyword'], 'fastq')

            del records[record[0]]

    output['with_keyword'].close()

    output['without_keyword'] = open(filename['without_keyword'], 'w')
    [SeqIO.write(record, output['without_keyword'], 'fastq')
     for (id, record) in records.items()]
    output['without_keyword'].close()

    print('\nCheck the results in ', filename[
          'without_keyword'], ' and ', filename['with_keyword'], '\n', sep='')

    end_time = time()

    print('\nTook %.3f seconds...\n' % (end_time - start_time))
