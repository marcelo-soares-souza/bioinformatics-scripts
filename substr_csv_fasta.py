#!/usr/bin/env python

import os
import sys
import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import Alphabet

if len(sys.argv) < 2:
    print('Usage: ', str(sys.argv[0]), '[CSV FILE] [FASTA FILE]')
else:

    records = SeqIO.index(str(sys.argv[2]), 'fasta')
    output = open(os.path.splitext(sys.argv[2])[0] + '.result.fasta', "w")

    with open(str(sys.argv[1]), 'r') as csv_file:
        reader = csv.reader(csv_file)

        for id, begin, end in reader:
            print('Processing', id, 'of size', len(records[id]), 'starting in', begin, 'ending', end)

            if begin < end:
                seq = records[id].seq[int(begin):int(end)]
            else:
                seq = records[id].seq[int(end):int(begin)].reverse_complement()

            header = '%s|size%s[%s_to_%s](%s nts)' % (id, len(records[id]), begin, end, len(seq))

            record = SeqRecord(Seq(str(seq), Alphabet()), id=str(header), description='')

            SeqIO.write(record, output, "fasta")

    print('Result in ', os.path.splitext(sys.argv[2])[0] + '.result.fasta')

    output.close()
    records.close()

