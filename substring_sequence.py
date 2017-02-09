#!/usr/bin/env python3

# (C) 2017 Marcelo Soares Souza <marcelo@riseup.net>
# This program is licensed under a LGPLv3 License.

import sys
import csv
import click

from time import time
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import Alphabet

@click.command()
@click.option('--csv', prompt='CSV File', help='CSV File')
@click.option('--input', prompt='FASTA or FASTQ File', help='FASTA or FASTQ File')
@click.option('--type', prompt='FASTA or FASTQ Type', help='--type FASTA or --type FASTQ')
@click.option('--output', prompt='Output File', help='Output File')
@click.option('--delimiter', prompt='Delimiter in CSV File', help='Delimiter in CSV File')

def substring_sequence(arg_csv, arg_input, arg_type, arg_output, arg_delimiter):
    start_time = time()

    csv_file = arg_csv
    records = SeqIO.index(arg_input, arg_type.lower())
    output = open(arg_output, "w")

    print('Using', csv_file, 'CSV Table and', arg_input, arg_type, ' File')

    with open(csv_file, 'r') as csv_file:
        reader = csv.reader(csv_file)
        reader = csv.reader(csv_file, delimiter=arg_delimiter)

        for id, begin, end in reader:

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

            header = '%s|size%s[%s_to_%s](%s nts)' % (
                id_name, len(records[id]), begin, end, len(seq))

            print(header)

            record = SeqRecord(Seq(str(seq), Alphabet()), id=str(header), description='')

            SeqIO.write(record, output, 'fasta')

    print('\nCheck the results in ', sys.argv[3], '\n', sep='')

    output.close()
    records.close()

    end_time = time()

    print('Took %.3f seconds...\n' % (end_time - start_time))

if __name__ == '__main__':
    substring_sequence()
