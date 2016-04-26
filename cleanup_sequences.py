#! /usr/bin/python3

# Cleanup Sequences Based on a .T6 File and a GI List (Optional)
# (C)2016 Marcelo Soares Souza <marcelo.soares@colaborador.embrapa.br>

import json
import re
from csv import reader
from datetime import datetime
from sys import argv
from os.path import splitext
from time import time
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import Alphabet

if len(argv) < 6:
    print('Usage:', str(argv[0]), '[T6 FILE] [INFO FILE] [FAST(A/Q) FILE] [PIDENT] [QCOVS] <GI LIST FILE>')
else:
    start_time = time()

    date = datetime.now().strftime('%Y%m%d')

    filename = {}
    filename['t6'] = str(argv[1])
    filename['info'] = str(argv[2])
    filename['data'] = str(argv[3])
    filename['data_format'] = str(splitext(filename['data'])[1][1:])
    filename['clean'] = str(splitext(filename['data'])[0] + '.clean.' + date + '.' + filename['data_format'])
    filename['removed'] = str(splitext(filename['data'])[0] + '.removed.' + date + '.' + filename['data_format'])

    params = {}
    params['pident'] = float(argv[4])
    params['qcovs'] = float(argv[5])

    config = json.loads(open(filename['info']).read())

    index = {}
    index['id'] = config['fields'].index('id')
    index['subject'] = config['fields'].index('subject')
    index['pident'] = config['fields'].index('pident')
    index['qcovs'] = config['fields'].index('qcovs')
    index['organism'] = config['fields'].index('organism')

    output = {}
    write = False
    use_gi = False
    removed_sequences = 0

    if len(argv) > 6:
        filename['gi'] = str(argv[6])
        use_gi = True

    print('\nReading', filename['data'], 'FAST(A/Q) File')
    records = SeqIO.to_dict(SeqIO.parse(filename['data'], filename['data_format']))

    print('\nUsing', filename['t6'], 'T6 Table and', filename['data'], 'FAST(A/Q) File')

    if use_gi:
        print('\nUsing', filename['gi'], 'GI List');

        removed_sequences_using_gi = 0

        filename['gi_output'] = str(splitext(filename['data'])[0] + '.gilist.' + date + '.' + filename['data_format'])
        output['gi_output'] = open(filename['gi_output'], 'w')

        with open(filename['gi'], 'r') as gi:
            gi_to_clean = {}
            gi_list = list(reader(gi))

            for gi in gi_list:
                gi_subject = 'gi|%s|' % (''.join(gi))
                gi_to_clean[gi_subject] = ''

    with open(filename['t6'], 'r') as t6:

        output['removed'] = open(filename['removed'], 'w')

        data = reader(t6, delimiter='\t')

        for value in data:
            id = value[index['id']]
            subject = value[index['subject']]
            pident = float(value[index['pident']])
            qcovs = float(value[index['qcovs']])
            organism = str(value[index['organism']])

            s = subject.rsplit('|',3)[0] + '|'

            if use_gi:
                if s in gi_to_clean.keys():
                    if id in records.keys():
                        print('\nFound GI Element', s, 'id', id)

                        header = '%s %s %s %s %s' % (id, subject, str(pident), str(qcovs), str(organism))
                        record = SeqRecord(Seq(str(records[id].seq), Alphabet()), id=str(header), description='')

                        if filename['data_format'] == 'fastq':
                            record.letter_annotations["phred_quality"] = records[id].letter_annotations["phred_quality"]

                        SeqIO.write(record, output['gi_output'], filename['data_format'])

                        del records[id]

                        removed_sequences_using_gi = removed_sequences_using_gi + 1

            if pident >= params['pident'] and qcovs >= params['qcovs']:
                if id in records.keys():

                    print('\nRemoving', id, '(', subject,')', organism, ', with PIDENT', pident, 'and QCOVS', qcovs)

                    header = '%s %s %s %s %s' % (id, subject, str(pident), str(qcovs), str(organism))
                    record = SeqRecord(Seq(str(records[id].seq), Alphabet()), id=str(header), description='')

                    if filename['data_format'] == 'fastq':
                        record.letter_annotations["phred_quality"] = records[id].letter_annotations["phred_quality"]

                    SeqIO.write(record, output['removed'], filename['data_format'])

                    del records[id]

                    removed_sequences = removed_sequences + 1

                    write = True

    output['removed'].close()

    if use_gi:
        output['gi_output'].close()

    print('\nWriting Results...')

    if write:
        output['clean'] = open(filename['clean'], 'w')
        [SeqIO.write(record, output['clean'], filename['data_format']) for (id, record) in records.items()]
        output['clean'].close()

        print('\nCheck the results in ', filename['clean'], '\n', sep='')
        print('Removed Sequences (', removed_sequences, ') in ', filename['removed'], '\n', sep='')

        if use_gi:
            print('Removed using GI List (', removed_sequences_using_gi, ') in ', filename['gi_output'], '\n', sep='')

    else:
        print('No sequences found\n')

    end_time = time()

    print('Took %.3f seconds...\n' % (end_time - start_time))
