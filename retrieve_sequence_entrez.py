#!/usr/bin/python3

# (C) 2016 Marcelo Soares Souza <marcelo.soares@colaborador.embrapa.br>
# This program is licensed under a LGPLv3 License.

import sys
from time import time
from Bio import Entrez
from Bio import SeqIO

Entrez.email = "lbb@cnpae.embrapa.br"

if len(sys.argv) < 2:
    print('Usage: ', str(sys.argv[0]), '[ID]')
else:
    start_time = time()

    id = str(sys.argv[1])

    handle = Entrez.efetch(db="nucleotide", id=id,
                           rettype="fasta", retmode="text")

    seq_record = SeqIO.read(handle, "fasta")

    output = open(id + '.fasta', 'w')

    SeqIO.write(seq_record, output, 'fasta')

    print("%s with %i features" % (seq_record.id, len(seq_record.features)))

    handle.close()
    output.close()

    end_time = time()

    print('Took %.3f seconds...\n' % (end_time - start_time))
