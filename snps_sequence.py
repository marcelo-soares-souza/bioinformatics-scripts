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

def verify_parameters(prefix):
    exit = 0

    if prefix == None:
        print('Please inform the Prefix --prefix option')
        exit = 1

    if exit == 1:
         sys.exit()

@click.command()
@click.option('--prefix', '-p', help='Prefix')

def snps_sequence(prefix):
    start_time = time()

    verify_parameters(prefix)

    files = ['bowtie_freebayes', 'bowtie_GATK', 'bowtie_samtools', 'bwa_freebayes', 'bwa_GATK', 'bwa_samtools']

    filename = {}
    output = {}

    for file in files:
        filename[file] = prefix + '_' + file

    for file in files:
        output[file] = open(filename[file] + '.csv', 'w')

    for key, value in filename.items():
        print("Processing", value)

        results = {}

        for line in open(value + '.vcf'):
            if line.startswith('#'):
                continue

            data = line.split();
            processed_line = ('%s\t%s\t%s\t%s\n' % (data[0] , data[1], data[3], data[4]))

            results[processed_line] = 1

        for result in results:
            output[key].write(result)

    print('Took %.3f seconds...\n' % (time() - start_time))

    for key, value in output.items():
        value.close()

if __name__ == '__main__':
    snps_sequence()


