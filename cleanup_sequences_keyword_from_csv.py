#!/usr/bin/python3

# (C) 2016 Marcelo Soares Souza <marcelo.soares@colaborador.embrapa.br>
# This program is licensed under a LGPLv3 License.

import os
import sys
import csv
import re
import operator
from time import time
from Bio import SeqIO

if len(sys.argv) < 3:
    print('Usage:', str(sys.argv[0]), '[FASTA FILE] [CSV FILE]')
else:
    start_time = time()

    fasta_filename = str(sys.argv[1])
    csv_filename = str(sys.argv[2])

    filename = {}
    filename['without_keyword'] = os.path.splitext(
        sys.argv[1])[0] + '.without' + '.fasta'
    filename['with_keyword'] = os.path.splitext(
        sys.argv[1])[0] + '.with' + '.fasta'

    print('\nReading', fasta_filename, 'FASTA File')

    records = SeqIO.to_dict(SeqIO.parse(fasta_filename, 'fasta'))

    output = {}
    output['with_keyword'] = open(filename['with_keyword'], 'w')

    records_sorted = sorted(records.items(), key=lambda v: int(v[0].split("_")[1]), reverse=False)

    with open(csv_filename, 'r') as csv_file:
        reader = csv.reader(csv_file)
        for keyword in reader:
            kw = ''.join(keyword) + '\n'

            for record in records_sorted:
                if kw.lower().find(record[1].description.lower()) != -1:
                    print('\nRemoving', record[1].description)
                    SeqIO.write(record[1], output['with_keyword'], 'fasta')
                    del records[record[0]]

    output['with_keyword'].close()

    output['without_keyword'] = open(filename['without_keyword'], 'w')
    [SeqIO.write(record, output['without_keyword'], 'fasta') for (id, record) in records_sorted]
    output['without_keyword'].close()

    print('\nCheck the results in ', filename[
          'without_keyword'], ' and ', filename['with_keyword'], '\n', sep='')

    end_time = time()

    print('\nTook %.3f seconds...\n' % (end_time - start_time))
