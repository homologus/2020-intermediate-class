#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Reads a file and applies all three frame shifts to
# it for all possible proteins a DNA sequence can
# code for. Programmed by Samuel Davenport

#from Bio.Seq import Seq
#from Bio.Alphabet import IUPAC
from Bio import SeqIO
from Bio.Blast import NCBIWWW
import os
from threading import Thread

NAME = 'Protien Shift'
__version__ = '0.0.1'

file = '/share/SARS/SARS-2020.fasta'
ftype = 'fasta'
wfile = 'output.txt'

progs = ['blastn', 'blastp', 'blastx', 'tblastn', 'tblastx']

if os.sys.argv[1:]:
    if '-f' in os.sys.argv:
        idx = os.sys.argv.index('-f')
        file = os.sys.argv[idx+1]
    if '-t' in os.sys.argv:
        idx = os.sys.argv.index('-t')
        ftype = os.sys.argv[idx+1]

def get_frames(dna_seq, l=100, n=3):
    """Returns protiens with a length at least <l> long from frame shifts <n>"""
    # For each frame shift,
    data = []
    for shift in range(abs(n) % 4):
        leng = len(dna_seq.seq[shift:])
        mod = leng % 3
        print(leng, mod)
        # Get the frame shift of the DNA sequence and
        # translate data to protiens
        protiens_data = dna_seq.seq[shift:-mod].translate()
        # Split the protiens by their end codons
        protiens = str(protiens_data).split('*')
        # For each protien we found,
##        for protien in protiens:
##            # If the protien has at least 100 amino acids,
##            if len(protien) >= l:
##                # Write the protien to the file
##                file.write(protien+'\n')
        data += protiens
    return [protien for protien in data if len(protien) >= l]

def search_blast(protien, numHits=50):
    result_handle = NCBIWWW.qblast("blastp", "nr", protien,
                                   hitlist_size=int(numHits), format_type='HTML')
    
    save_file = open("my_blast.xml", "w")
    data = result_handle.read()
    #text = data.split('<Iteration>')[1].split('</Iteration_hits>')[0]
    #text = ' '.join([i for i in ' '.join([i for i in text.split('\n')]).split(' ') if i != ''])
    
    save_file.write(data)#result_handle.read())
    save_file.close()
    result_handle.close()
    
    text = [i.split('</Hit_def>\n')[0] for i in data.split('</Hit_id>\n')][1:]
    names = [i.split('  <Hit_def>')[1] for i in text]
    return names

class blast_thread(Thread):
    def __init__(self, protien, dataList, numHits=5):
        Thread.__init__(self)
        self.protien = protien
        self.dataList = dataList
        self.hits = numHits
        self.start()
    
    def run(self):
        self.dataList.append(search_blast(self.protien, self.hits))
    pass

def get_data_from_blast(protiens, numHitsEach=5):
    data = []
    threads = {}
    for i in range(len(protiens)):
        protien = protiens[i]
        threads[i] = blast_thread(protien, data, numHitsEach)

def run():
    print('Opening file', file, 'type', ftype)
    
    dna_seq = SeqIO.read(file, ftype)
##    dna_seqs = SeqIO.parse(file, ftype)
    file = open(wfile, 'w')
##    for dna_seq in dna_seqs:
    protiens = get_frames(dna_seq, l=100, n=3)
    file.write('\n'.join(protiens)+'\n\n')
    names = ['\n'.join(protiens)+'\n\n' for names in [search_blast(prot, 5) for prot in protiens]]
    file.write(''.join([str(i)+': '+names[i] for i in range(len(names))]))
    print('Done')

if __name__ != '__main__':
    run()

