#!/usr/bin/env python3

# (C) 2017 Marcelo Soares Souza <marcelo@riseup.net>
# This program is licensed under a LGPLv3 License.

import sys
import csv
import click

from time import time
from collections import defaultdict

from Bio import SeqIO


def verify_parameters(prefix):
    exit = 0

    if prefix is None:
        print('Please inform the Prefix --prefix option')
        exit = 1

    if exit == 1:
        sys.exit()


def clean_vcfs(prefix, files, filenames):
    output = {}

    for file in files:
        output[file] = open(filenames[file] + '.csv', 'w')

    for key, value in filenames.items():
        print('Cleaning VCF:', value)

        results = {}

        for line in open(value + '.vcf'):
            if line.startswith('#'):
                continue

            data = line.split()
            processed_key = ('%s\t%s' % (data[0], data[1]))
            processed_val = ('%s\t%s\t%s\t%s\n' % (data[0], data[1], data[3], data[4]))

            results[processed_key] = processed_val

        for k, v in results.items():
            output[key].write(results[k])

    for key, value in output.items():
        value.close()


def process(prefix, filenames):
    sequences = defaultdict(lambda: defaultdict(list))

    for suffix, filename in filenames.items():

        with open(filename + '.csv', 'r') as csv_data:
            reader = csv.reader(csv_data, delimiter='\t')

            for data in reader:
                sequence = data[0]
                position = data[1]

                sequences[sequence][position].append(suffix)

    tools = defaultdict(list)

    for sequence, lists in sequences.items():
        for key, suffix in lists.items():
            tools[', '.join(suffix)].append(sequence + ' ' + key)

    tools_sorted = sorted(tools, key=len, reverse=True)

    output = open(prefix + '.result', 'w')

    for tool in tools_sorted:
        header = ','.join(sorted(tool.strip().split(','), key=len, reverse=True)).strip()

        output.write('Found on %s - %s\n\n' % (header, len(tools[tool])))

        for seq_pos in tools[tool]:
            output.write('%s\n' % (seq_pos))

        output.write('\n')

    print('\nIntermediate result saved in', prefix + '.result\n')

    output.close()


def process_fasta(fasta, updown, prefix, filenames):
    organism = prefix.split('_')[0]

    output = {}
    fasta_file = fasta
    up_down_value = updown

    chrs = defaultdict(lambda: defaultdict(list))
    sequences = {}

    records = SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta'))

    for suffix, filename in filenames.items():

        with open(filename + '.csv', 'r') as csv_data:
            reader = csv.reader(csv_data, delimiter='\t')

            for data in reader:
                chr = data[0]
                pos = data[1]
                ref = data[2]
                alt = data[3]
                snp_id = organism + '_' + chr + '_' + pos

                begin = int(pos) - up_down_value
                end = int(pos) + up_down_value

                if begin > 0:
                    seq = str(records[chr].seq[int(begin):int(end)])
                    seq = seq[:up_down_value] + '[' + ref + '/' + alt + ']' + seq[up_down_value:]
                    info = '%s;%s;%s;%s;%s;%s' % (snp_id, seq, chr, pos, ref, alt)
                    seq_pos = '%s %s' % (chr, pos)
                    sequences[seq_pos] = info
                    chrs[chr][pos].append(suffix)

    tools = defaultdict(list)

    for chr, lists in chrs.items():
        for key, suffix in lists.items():
            tools[', '.join(suffix)].append(chr + ' ' + key)

    tools_sorted = sorted(tools, key=len, reverse=True)

    seqs = {}

    for tool in tools_sorted:
        len_tool = str(len(tool.strip().split(',')))

        for seq_pos in tools[tool]:
            seqs[seq_pos] = '%s;%s' % (len_tool, (sequences[seq_pos]))

    last_pos = pos = 0
    last_chr = chr = ''

    for i in range(1, 7):
        output[i] = open(prefix + '.' + str(i) + '.csv', 'w')

    output[0] = open(prefix + '.csv', 'w')

    for k in sorted(seqs, key=lambda x: (x.split()[0], int(x.split()[1]))):

        tool = int(seqs[k].split(';')[0])
        chr = str(seqs[k].split(';')[3])
        pos = int(seqs[k].split(';')[4])

        col = '%s\n' % (seqs[k])

        if last_pos != 0 and last_chr == chr:
            col = '%s;%s\n' % (seqs[k], str((pos - last_pos)))

        output[tool].write(col)
        output[0].write(col)

        last_chr = chr
        last_pos = pos

    for i in range(0, 7):
        output[i].close()


@click.command()
@click.option('--prefix', '-p', help='Prefix')
@click.option('--fasta', '-f',  help='FASTA File')
@click.option('--updown', '-u', help='Up/Down Stream Value', default=35)
def snps_sequence(prefix, fasta, updown):
    start_time = time()

    verify_parameters(prefix)

    files = ['bowtie_freebayes', 'bowtie_GATK', 'bowtie_samtools', 'bwa_freebayes', 'bwa_GATK', 'bwa_samtools']
    filenames = {}

    for file in files:
        filenames[file] = prefix + '_' + file

    clean_vcfs(prefix, files, filenames)

    process(prefix, filenames)

    if fasta is not None:
        process_fasta(fasta, updown, prefix, filenames)

    print('Took %.3f seconds...\n' % (time() - start_time))


if __name__ == '__main__':
    snps_sequence()
