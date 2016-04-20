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

if len(sys.argv) < 3:
    print('Usage:', str(sys.argv[0]), '[CSV FILE] [PIDENT]')
else:
    start_time = time()

    csv_filename = str(sys.argv[1])
    arg_pident = float(sys.argv[2])

    print('Using', csv_filename, '\n')

    result = defaultdict(list)
    info = defaultdict(dict)

    with open(csv_filename, 'r') as csv_data:
        reader = csv.reader(csv_data, delimiter='\t')

        for data in reader:
            pident = float(data[11])

            if pident >= arg_pident:
                id = str(data[2])
                slen = data[3]
                start, end  = int(data[6]), int(data[7])

                if start > end:
                    start, end = end, start

                # print("id", id, "slen", slen, "start", start, "end", end, "pident", pident)
                info[id]['slen'] = slen
                info[id]['pident'] = pident

                result[id].append([start, end])

    result = sorted(dict(result).items())

    output_filename = os.path.splitext(sys.argv[1])[0] + '.result-PIDENT_' + str(arg_pident) + '.' + datetime.datetime.now().strftime('%Y%m%d') + '.txt'
    output = open(output_filename, 'w')

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

        log = 'SSEQID: %s SLEN: %s BPs: %s PIDENT >= %s\n' % (str(k), str(info[k]['slen']), str(size), str(arg_pident))
        print(log)
        output.write(log)

    output.close()

    end_time = time()

    print('\nCheck the results in ', output_filename, '\n', sep='')
    print('\nTook %.3f seconds...\n' % (end_time - start_time))
