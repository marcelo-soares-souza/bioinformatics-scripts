#!/usr/bin/python3

# (C) 2016 Marcelo Soares Souza <marcelo.soares@colaborador.embrapa.br>
# This program is licensed under a LGPLv3 License.

import os
import sys
import pandas
from concurrent.futures import ProcessPoolExecutor
from time import time
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import Alphabet


def slicing(data):
    id, begin, end, record = data['id'], data[
        'begin'], data['end'], data['record']

    print('\nProcessing', id, 'of size', len(record),
          'starting in', begin, 'ending', end)

    if begin < end:
        print('Normal Slicing...')
        seq = record.seq[begin:end]
    else:
        print('Reverse Complement...')
        seq = record.seq[end:begin].reverse_complement()

    header = '%s|size%s[%s_to_%s](%s nts)' % (
        id, len(record), begin, end, len(seq))

    return SeqRecord(Seq(str(seq), Alphabet()), id=str(header), description='')


def main():
    if len(sys.argv) < 2:
        print('Usage: ', str(sys.argv[0]), '[CSV FILE] [FASTA FILE]')
    else:
        start_time = time()

        csv_file = str(sys.argv[1])
        fasta_file = str(sys.argv[2])

        print('\nUsing CSV', csv_file, 'and FASTA', fasta_file)

        records = SeqIO.index(fasta_file, 'fasta')
        output = open(os.path.splitext(sys.argv[2])[0] + '.result.fasta', 'w')
        pool = ProcessPoolExecutor()
        results = []
        data = []

        print('Reading CSV and FASTA File...')

        create_start_time = time()

        reader = list(pandas.read_csv(
            csv_file, header=None, index_col=None).values)

        for id, begin, end in reader:
            register = {}
            register['id'] = id
            register['begin'] = int(begin)
            register['end'] = int(end)
            register['record'] = records[id]

            data.append(dict(register))

        create_end_time = time()

        print('Create Registers Time:', (create_end_time - create_start_time))

        print('Processing...')

        results = list(pool.map(slicing, data))

        [SeqIO.write(result, output, 'fasta') for result in results]

        print('\nCheck the results in ', os.path.splitext(
            sys.argv[2])[0], '.result.fasta\n', sep='')

        output.close()
        records.close()

        end_time = time()

        print('Took %.3f seconds...\n' % (end_time - start_time))

if __name__ == '__main__':
    main()
