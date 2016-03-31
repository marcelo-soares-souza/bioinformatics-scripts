#!/usr/bin/python3

import os
import sys
import csv
import concurrent.futures
from concurrent.futures import ProcessPoolExecutor
from time import time
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import Alphabet


def slicing(record):
    print('\nProcessing', record['id'], 'of size', len(record['record']), 'starting in', record['begin'], 'ending', record['end'])

    if record['begin'] < record['end']:
        print('Normal Slicing...')
        seq = record['record'].seq[record['begin']:record['end']]
    else:
        print('Reverse Complement...')
        seq = record['record'].seq[record['end']:record['begin']].reverse_complement()

    header = '%s|size%s[%s_to_%s](%s nts)' % (record['id'], len(record['record']), record['begin'], record['end'], len(seq))

    return SeqRecord(Seq(str(seq), Alphabet()), id=str(header), description='')

if len(sys.argv) < 2:
    print('Usage: ', str(sys.argv[0]), '[CSV FILE] [FASTA FILE]')
else:
    start_time = time()

    csv_file = str(sys.argv[1])
    fasta_file = str(sys.argv[2])

    records = SeqIO.index(fasta_file, 'fasta')
    output = open(os.path.splitext(sys.argv[2])[0] + '.result.fasta', 'w')
    results = []
    datas = []

    pool = ProcessPoolExecutor(max_workers=1)

    print('Using', csv_file, 'CSV Table and', fasta_file, 'FASTA File')


    with open(csv_file, 'r') as csv_file:
        reader = csv.reader(csv_file)

        for id, begin, end in reader:
            data = {}
            data['id'] = id
            data['record'] = records[id]
            data['begin'] = int(begin)
            data['end'] = int(end)

            datas.append(dict(data))

    results = list(pool.map(slicing, datas))

    [SeqIO.write(result, output, 'fasta') for result in results]

    print('\nCheck the results in ', os.path.splitext(sys.argv[2])[0], '.result.fasta\n', sep='')

    output.close()
    records.close()

    end_time = time()

    print('Took %.3f seconds...\n' % (end_time - start_time))
