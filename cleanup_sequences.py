#! /usr/bin/python3

# Cleanup Sequences Based on a .T6 File and a GI List (Optional)
# (C)2016 Marcelo Soares Souza <marcelo.soares@colaborador.embrapa.br>

from json import loads
from yaml import load, dump
from csv import reader
from datetime import datetime
from sys import argv
from os.path import splitext
from time import time
from collections import defaultdict
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

    try:
        config = loads(open(filename['info']).read())
    except:
        print('\nTrying YAML Based File')
        config = load(open(filename['info']).read())
    else:
        print('\nLoad JSON Based File')

    params = {}
    params['pident'] = float(config['pident'])
    params['qcovs'] = int(config['qcovs'])

    filename['t6'] = str(config['t6'])
    filename['data'] = str(config['input'])
    filename['data_format'] = str(splitext(filename['data'])[1][1:])
    filename['clean'] = str(splitext(filename['data'])[0] +
                            '.nohits-' +
                            str(config['blast-type']).lower() +
                            '-' +
                            str(config['blast-db']).lower() +
                            '.' + filename['data_format'])

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
    removed_extension = ''

    if 'gi' in config:
        filename['gi'] = str(config['gi'])
        params['gi-type'] = str(config['gi-type']).lower()

        if params['gi-type'] == 'negative':
          removed_extension = '_neg-list-%s' % str(config['gi-taxon']).lower()
        else:
          removed_extension = '_pos-list-%s' % str(config['gi-taxon']).lower()

        use_gi = True

    filename['removed'] = str(splitext(filename['data'])[0] +
                              '.filter_pident-' +
                              str(params['pident'])  +
                              '_qcovs-' +
                              str(params['qcovs']) +
                              removed_extension + 
                              '.' + filename['data_format'])

    filename['removed-stats'] = str(splitext(filename['removed'])[0]) + '.stats'

    print('\nReading', filename['data'], 'FAST(A/Q) File')
    records = SeqIO.to_dict(SeqIO.parse(filename['data'], filename['data_format']))

    print('\nUsing', filename['t6'], 'T6 Table and', filename['data'], 'FAST(A/Q) File')

    if use_gi:
        print('\nUsing', filename['gi'], 'GI List as', params['gi-type'].title());

        removed_sequences_using_gi = 0

        gi_output_extension = ''

        if params['gi-type'] == 'negative':
            gi_output_extension = '.at-gilist-'
        else:
            gi_output_extension = '.not-gilist-'

        filename['gi_output'] = str(splitext(filename['data'])[0] +
                                    gi_output_extension +
                                    str(config['gi-taxon']).lower() +
                                    '.' + filename['data_format'])

        output['gi_output'] = open(filename['gi_output'], 'w')

        with open(filename['gi'], 'r') as gi:
            gi_to_clean = {}
            gi_list = list(reader(gi))

            for gi in gi_list:
                gi_sseqid = 'gi|%s|' % (''.join(gi))
                gi_to_clean[gi_sseqid] = ''

    stats = defaultdict(int)

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
                if s in gi_to_clean.keys() and params['gi-type'] == 'negative':
                    if qseqid in records.keys():
                        print('\nFound GI Element', s, 'qseqid', qseqid, sscinames)

                        header = '%s %s %s %s %s' % (qseqid, sseqid, str(pident), str(qcovs), str(sscinames))
                        record = SeqRecord(Seq(str(records[qseqid].seq), Alphabet()), id=str(header), description='')

                        if filename['data_format'] == 'fastq':
                            record.letter_annotations["phred_quality"] = records[qseqid].letter_annotations["phred_quality"]

                        SeqIO.write(record, output['gi_output'], filename['data_format'])

                        del records[qseqid]

                        removed_sequences_using_gi = removed_sequences_using_gi + 1

                        write = True


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

                    stats[sscinames] = stats[sscinames] + 1

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
        print('Stats for Removed Sequences in ', filename['removed-stats'], '\n', sep='')

        if use_gi:
            print('Removed using GI List (', removed_sequences_using_gi, ') in ', filename['gi_output'], '\n', sep='')

    else:
        print('No sequences found\n')

    end_time = time()

    print('Took %.3f seconds...\n' % (end_time - start_time))

    total_stats = sum(stats.values())
    output['removed-stats'] = open(filename['removed-stats'], 'w')
    [output['removed-stats'].write("%s: %s of %s (%.2f%%) \n" % (k, v, total_stats, float(v * 100 / total_stats))) for k, v in sorted(stats.items(), key=lambda x: x[1], reverse=True)]
    output['removed-stats'].close()
