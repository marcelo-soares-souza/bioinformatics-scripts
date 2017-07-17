#!/usr/bin/env python3

# (C) 2017 Marcelo Soares Souza <marcelo@riseup.net>
# This program is licensed under a LGPLv3 License.

import sys
from csv import reader
from time import time

import click
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
    output_io = open(output, "w")

    print('Using', csv, 'CSV Table and', input, format.upper(), ' File')

    with open(csv, 'r') as csv_io:
        data = reader(csv_io, delimiter=';')

        for id in data:

            id = "".join(id)

            print('\nProcessing', id)

            # seq = records[id]

            # print(str(seq))

            # record = SeqRecord(Seq(str(seq), Alphabet()), id=str(id), description='')
            record = records[id]
            SeqIO.write(record, output_io, 'fasta')

    print('\nCheck the results in ', output, '\n', sep='')

    output_io.close()
    records.close()

    print('Took %.3f seconds...\n' % (time() - start_time))


if __name__ == '__main__':
    substring_sequence()
