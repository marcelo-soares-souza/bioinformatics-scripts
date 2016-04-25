#! /usr/bin/python3

# Cleanup Sequences Based on a .T6 File and a GI List (Optional)
# (C)2016 Marcelo Soares Souza <marcelo.soares@colaborador.embrapa.br>

from csv import reader
from datetime import datetime
from sys import argv
from os.path import splitext
from time import time
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import Alphabet

if len(argv) < 5:
    print('Usage:', str(argv[0]), '[T6 FILE] [FAST(A/Q) FILE] [PIDENT] [QCOVS] [GI LIST FILE]')
else:
    start_time = time()

    date = datetime.now().strftime('%Y%m%d')

    filename = {}
    filename['t6'] = str(argv[1])
    filename['data'] = str(argv[2])
    filename['data_format'] = str(splitext(filename['data'])[1][1:])
    filename['clean'] = str(splitext(filename['data'])[0] + '.clean.' + date + '.' + filename['data_format'])
    filename['removed'] = str(splitext(filename['data'])[0] + '.removed.' + date + '.' + filename['data_format'])

    if len(argv) > 6:
        filename['gi'] = str(argv[5])
	
    params = {}
    params['pident'] = float(argv[3])
    params['qcovs'] = float(argv[4])

    output = {}
    write = False
    removed_sequences = 0

    print('\nReading', filename['data'], 'FAST(A/Q) File')
    records = SeqIO.to_dict(SeqIO.parse(filename['data'], filename['data_format']))

    print('\nUsing', filename['t6'], 'T6 Table and', filename['data'], 'FAST(A/Q) File')

    with open(filename['t6'], 'r') as t6:
        data = reader(t6, delimiter='\t')

        output['removed'] = open(filename['removed'], 'w')

        for value in data:
            id = value[0]
            subject = value[2]
            pident = float(value[11])
            qcovs = float(value[14])

            if pident >= params['pident'] and qcovs >= params['qcovs']:
                if id in records.keys():

                    print('\nRemoving', id, 'with PIDENT', pident, 'and QCOVS', qcovs)

                    header = '%s %s %s %s' % (id, subject, str(pident), str(qcovs))
                    record = SeqRecord(Seq(str(records[id].seq), Alphabet()), id=str(header), description='')

                    if filename['data_format'] == 'fastq':
                        record.letter_annotations["phred_quality"] = records[id].letter_annotations["phred_quality"]

                    SeqIO.write(record, output['removed'], filename['data_format'])

                    del records[id]

                    write = True
                    removed_sequences = removed_sequences + 1

    output['removed'].close()

    print('\nWriting Results...')

    if write:
        output['clean'] = open(filename['clean'], 'w')
        [SeqIO.write(record, output['clean'], filename['data_format']) for (id, record) in records.items()]
        output['clean'].close()

        print('\nCheck the results in ', filename['clean'], '\n', sep='')
        print('Removed Sequences (', removed_sequences, ') in ', filename['removed'], '\n', sep='')
    else:
        print('No sequences found\n')

    end_time = time()

    print('Took %.3f seconds...\n' % (end_time - start_time))
