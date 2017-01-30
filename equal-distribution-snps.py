#!/usr/bin/env python3

# (C) 2016 Marcelo Soares Souza <marcelo@riseup.net>
# This program is licensed under a LGPLv3 License.

import sys
import csv
from collections import defaultdict
from time import time
from math import ceil


def takespread(sequence, num):
    length = float(len(sequence))

    for i in range(num):
        yield sequence[int(ceil(i * length / num))]


if len(sys.argv) < 2:
    print('Usage:', str(sys.argv[0]), '[CSV FILE]')
else:
    start_time = time()

    csv_filename = str(sys.argv[1])

    print('Using', csv_filename, '\n')

    info = defaultdict(list)
    snp = defaultdict(dict)

    with open(csv_filename, 'r') as csv_data:
        reader = csv.reader(csv_data, delimiter=',')

        for data in reader:
            chr = data[1]
            pos = int(data[2])
            snp[pos] = data[0]

            info[chr].append(pos)

    info = sorted(dict(info).items())

    for k, v in info:
        print('Chromosome:', k, '\n')
        data = takespread(v, 10)

        for i in data:
            print("SNP: %s" % (snp[i]))

        print('\n')
