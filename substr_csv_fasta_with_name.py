#!/usr/bin/python3

# (C) 2016 Marcelo Soares Souza <marcelo.soares@colaborador.embrapa.br>
# This program is licensed under a LGPLv3 License.

import sys
import csv
from time import time
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import Alphabet

if len(sys.argv) < 4:
    print('Usage:', str(sys.argv[0]), '[CSV FILE] [FASTA FILE] [OUTPUT]')
else:
    start_time = time()

    csv_file = str(sys.argv[1])
    fasta_file = str(sys.argv[2])

    records = SeqIO.index(fasta_file, 'fasta')
    output = open(sys.argv[3], "w")

    print('Using', csv_file, 'CSV Table and', fasta_file, 'FASTA File')

    with open(csv_file, 'r') as csv_file:
        reader = csv.reader(csv_file, delimiter=';')

        for name, id, begin, end in reader:

            position_begin = int(begin) - 1

            print('\nProcessing', id, 'of size', len(
                records[id]), 'starting in', begin, 'ending', end)

            if int(begin) < int(end):
                print('Normal Slicing...')
                seq = records[id].seq[int(position_begin):int(end)]
                id_name = id
            else:
                print('Reverse Complement...')
                seq = records[id].seq[int(end):int(position_begin)].reverse_complement()
                id_name = id + '_RC'

            # header = '%s|size%s[%s_to_%s](%s nts)' % (
            #    id_name, len(records[id]), begin, end, len(seq))

            header = '%s' % (name)

            print(header)

            record = SeqRecord(Seq(str(seq), Alphabet()),
                               id=str(header), description='')

            SeqIO.write(record, output, 'fasta')

    print('\nCheck the results in ', sys.argv[3], '\n', sep='')

    output.close()
    records.close()

    end_time = time()

    print('Took %.3f seconds...\n' % (end_time - start_time))
