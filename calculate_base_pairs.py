#!/usr/bin/python3

import os
import sys
import csv
import datetime
from collections import defaultdict
from time import time
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import Alphabet

if len(sys.argv) < 2:
    print('Usage:', str(sys.argv[0]), '[CSV FILE]')
else:
    start_time = time()

    csv_filename = str(sys.argv[1])

    print('Using', csv_filename)

    result = defaultdict(list)

    with open(csv_filename, 'r') as csv_data:
        reader = csv.reader(csv_data, delimiter='\t')

        for data in reader:
            id = str(data[2])
            start, end  = int(data[5]), int(data[6])

            if start > end:
              start, end = end, start

            result[id].append([start, end])

    result = sorted(dict(result).items())

    for k, v in result:
        s = sorted(v)

        min_value = min(sum(s, []))
        max_value = max(sum(s, []))

        x = set(range(min_value, max_value))

        initial_size = len(x)

        for y in s:
          y = set(range(y[0], y[1]))
          x = x - y

        size = initial_size - len(x)

        print(k, size)

    end_time = time()

    print('\nTook %.3f seconds...\n' % (end_time - start_time))
