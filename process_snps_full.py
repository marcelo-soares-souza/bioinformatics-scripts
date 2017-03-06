#!/usr/bin/env python3

# (C) 2016 Marcelo Soares Souza <marcelo@riseup.net>
# This program is licensed under a LGPLv3 License.

import sys
import csv
from collections import defaultdict
from time import time
from Bio import SeqIO


def getKey(item):
    return item[1]


if len(sys.argv) < 4:
    print('Usage:', str(sys.argv[0]), '[PREFIX FILENAME] [FASTA] [UP/DO VALUE]')
else:
    start_time = time()

    prefix = str(sys.argv[1])

    filename = {}
    filename['bowtie_freebayes'] = prefix + '_' + 'bowtie_freebayes.csv'
    filename['bowtie_GATK'] = prefix + '_' + 'bowtie_GATK.csv'
    filename['bowtie_samtools'] = prefix + '_' + 'bowtie_samtools.csv'
    filename['bwa_freebayes'] = prefix + '_' + 'bwa_freebayes.csv'
    filename['bwa_GATK'] = prefix + '_' + 'bwa_GATK.csv'
    filename['bwa_samtools'] = prefix + '_' + 'bwa_samtools.csv'
    organism = prefix.split('_')[0]

    output = {}

    fasta_file = str(sys.argv[2])
    up_down_value = int(sys.argv[3])

    chrs = defaultdict(lambda: defaultdict(list))
    sequences = {}

    print('Using Prefix', prefix, '\n')

    records = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

    for suffix, filename in filename.items():

        with open(filename, 'r') as csv_data:
            reader = csv.reader(csv_data, delimiter='\t')

            for data in reader:
                chr = data[0]
                pos = data[1]
                ref = data[2]
                alt = data[3]
                snp_id = organism + "_" + chr + "_" + pos

                begin = int(pos) - up_down_value
                end = int(pos) + up_down_value

                if begin > 0:
                    seq = str(records[chr].seq[int(begin):int(end)])
                    seq = seq[:up_down_value] + '[' + ref + '/' + alt + ']' + seq[up_down_value:]
                    info = '%s;%s;%s;%s;%s;%s' % (snp_id, seq, chr, pos, ref, alt)
                    seq_pos = "%s %s" % (chr, pos)
                    sequences[seq_pos] = info
                    chrs[chr][pos].append(suffix)

    tools = defaultdict(list)

    for chr, list in chrs.items():
        for key, suffix in list.items():
            tools[', '.join(suffix)].append(chr + ' ' + key)

    tools_sorted = sorted(tools, key=len, reverse=True)

    seqs = {}

    for tool in tools_sorted:
        len_tool = str(len(tool.strip().split(',')))

        for seq_pos in tools[tool]:
            seqs[seq_pos] = "%s;%s" % (len_tool, (sequences[seq_pos]))

    last_pos = pos = 0
    last_chr = chr = ''

    output[6] = open(prefix + '.6.csv', 'w')
    output[5] = open(prefix + '.5.csv', 'w')
    output[4] = open(prefix + '.4.csv', 'w')
    output[3] = open(prefix + '.3.csv', 'w')
    output[2] = open(prefix + '.2.csv', 'w')
    output[1] = open(prefix + '.1.csv', 'w')
    output[0] = open(prefix + '.csv', 'w')

    for k in sorted(seqs, key=lambda x: (x.split()[0], int(x.split()[1]))):

        tool = int(seqs[k].split(';')[0])
        chr = str(seqs[k].split(';')[3])
        pos = int(seqs[k].split(';')[4])

        col = "%s\n" % (seqs[k])

        if last_pos != 0 and last_chr == chr:
            col = "%s;%s\n" % (seqs[k], str((pos - last_pos)))

        output[tool].write(col)
        output[0].write(col)

        last_chr = chr
        last_pos = pos

    output[6].close()
    output[5].close()
    output[4].close()
    output[3].close()
    output[2].close()
    output[1].close()

    end_time = time()

    print('\nTook %.3f seconds...\n' % (end_time - start_time))
