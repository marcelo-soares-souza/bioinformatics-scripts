#! /usr/bin/python3

# Cleanup Sequences Based on a .T6 File and a GI List (Optional)
# (C)2016 Marcelo Soares Souza <marcelo.soares@colaborador.embrapa.br>

# from json import loads
from yaml import load, dump
from csv import reader
from datetime import datetime
from sys import argv
from os.path import splitext
from time import time
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import Alphabet

if len(argv) < 2:
    print('Usage:', str(argv[0]), '[INFO FILE]')
else:
    start_time = time()

    date = datetime.now().strftime('%Y%m%d')

    filename = {}
    filename['info'] = str(argv[1])
    config = load(open(filename['info']).read())
    print(config)
    filename['t6'] = str(config['t6'])
    filename['data'] = str(config['input'])
    filename['data_format'] = str(splitext(filename['data'])[1][1:])
    filename['clean'] = str(splitext(filename['data'])[0] + '.clean.' + date + '.' + filename['data_format'])
    filename['removed'] = str(splitext(filename['data'])[0] + '.removed.' + date + '.' + filename['data_format'])


    params = {}
    params['pident'] = float(config['pident'])
    params['qcovs'] = int(config['qcovs'])

    index = {}
    index['qseqid'] = config['fields'].split().index('qseqid')
    index['sseqid'] = config['fields'].split().index('sseqid')
    index['pident'] = config['fields'].split().index('pident')
    index['qcovs'] = config['fields'].split().index('qcovs')
    index['sscinames'] = config['fields'].split().index('sscinames')

    output = {}
    write = False
    use_gi = False
    removed_sequences = 0

    if 'gi' in config:
        filename['gi'] = str(config['gi'])
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
                gi_sseqid = 'gi|%s|' % (''.join(gi))
                gi_to_clean[gi_sseqid] = ''

    with open(filename['t6'], 'r') as t6:

        output['removed'] = open(filename['removed'], 'w')

        data = reader(t6, delimiter='\t')

        for value in data:
            qseqid = value[index['qseqid']]
            sseqid = value[index['sseqid']]
            pident = float(value[index['pident']])
            qcovs = int(value[index['qcovs']])
            sscinames = str(value[index['sscinames']])

            s = sseqid.rsplit('|',3)[0] + '|'

            if use_gi:
                if s in gi_to_clean.keys():
                    if qseqid in records.keys():
                        print('\nFound GI Element', s, 'qseqid', qseqid)

                        header = '%s %s %s %s %s' % (qseqid, sseqid, str(pident), str(qcovs), str(sscinames))
                        record = SeqRecord(Seq(str(records[qseqid].seq), Alphabet()), id=str(header), description='')

                        if filename['data_format'] == 'fastq':
                            record.letter_annotations["phred_quality"] = records[qseqid].letter_annotations["phred_quality"]

                        SeqIO.write(record, output['gi_output'], filename['data_format'])

                        del records[qseqid]

                        removed_sequences_using_gi = removed_sequences_using_gi + 1

            if pident >= params['pident'] and qcovs >= params['qcovs']:
                if qseqid in records.keys():

                    print('\nRemoving', qseqid, '(', sseqid,')', sscinames, ', with PIDENT', pident, 'and QCOVS', qcovs)

                    header = '%s %s %s %s %s' % (qseqid, sseqid, str(pident), str(qcovs), str(sscinames))
                    record = SeqRecord(Seq(str(records[qseqid].seq), Alphabet()), id=str(header), description='')

                    if filename['data_format'] == 'fastq':
                        record.letter_annotations["phred_quality"] = records[qseqid].letter_annotations["phred_quality"]

                    SeqIO.write(record, output['removed'], filename['data_format'])

                    del records[qseqid]

                    removed_sequences = removed_sequences + 1

                    write = True

    output['removed'].close()

    if use_gi:
        output['gi_output'].close()

    print('\nWriting Results...')

    if write:
        output['clean'] = open(filename['clean'], 'w')
        [SeqIO.write(record, output['clean'], filename['data_format']) for (qseqid, record) in records.items()]
        output['clean'].close()

        print('\nCheck the results in ', filename['clean'], '\n', sep='')
        print('Removed Sequences (', removed_sequences, ') in ', filename['removed'], '\n', sep='')

        if use_gi:
            print('Removed using GI List (', removed_sequences_using_gi, ') in ', filename['gi_output'], '\n', sep='')

    else:
        print('No sequences found\n')

    end_time = time()

    print('Took %.3f seconds...\n' % (end_time - start_time))
