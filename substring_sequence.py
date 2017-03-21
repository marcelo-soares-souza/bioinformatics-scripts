#!/usr/bin/env python3

# (C) 2017 Marcelo Soares Souza <marcelo@riseup.net>
# This program is licensed under a LGPLv3 License.

import sys
# import csv
import click

from time import time
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import Alphabet


def verify_parameters(csv, input, output):
    exit = 0

    if csv is None:
        print("Please inform the CSV File with --csv option")
        exit = 1

    if input is None:
        print("Please inform the Input File with --input option")
        exit = 1

    if output is None:
        print("Please inform the Output File with --output option")
        exit = 1

    if exit == 1:
        sys.exit()


@click.command()
@click.option('--csv', '-c', help='CSV File')
@click.option('--input', '-i', help='FASTA/FASTQ File')
@click.option('--output', '-o', help='Output File')
@click.option('--format', '-f', default='fasta', help='--format fasta or --format fastq')
@click.option('--delimiter', '-d', default=';', help='Delimiter in CSV File')
def substring_sequence(csv, input, output, format,  delimiter):
    verify_parameters(csv, input, output)

    start_time = time()

    records = SeqIO.index(input, format.lower())
    output = open(output, "w")

    print('Using', csv, 'CSV Table and', input, format, ' File')

    with open(csv, 'r') as csv_file:
        reader = csv.reader(csv_file)
        reader = csv.reader(csv_file, delimiter=delimiter)

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

    print('Took %.3f seconds...\n' % (time() - start_time))


if __name__ == '__main__':
    substring_sequence()
