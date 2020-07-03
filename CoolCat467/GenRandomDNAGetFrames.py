#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Reads a file and applies all three frame shifts to
# it for all possible proteins a DNA sequence can
# code for. Programmed by Samuel Davenport

from Bio.Seq import Seq
from Bio import SeqIO
from random import randint, choice
import os
from threading import Thread

NAME = 'Threaded Random DNA and Shifts to Proteins'
__version__ = '0.0.2'

def get_frames(dna_seq, n=3):
    """Returns proteins with a length at least <l> long from frame shifts <n>"""
    # For each frame shift,
    data = []
    for shift in range(abs(n) % 4):
        leng = len(dna_seq[shift:])
        mod = leng % 3
        #print(leng, mod)
        # Get the frame shift of the DNA sequence and
        # translate data to proteins
        proteins_data = dna_seq[shift:-mod].translate()
        # Split the proteins by their end codons
        proteins = str(proteins_data).split('*')
        # For each protein we found,
##        for protein in proteins:
##            # If the protein has at least 100 amino acids,
##            if len(protein) >= l:
##                # Write the protein to the file
##                file.write(protein+'\n')
        data += proteins
    return [protein for protein in data if len(protein) >= 5]

def gen_random_proper_seq(length, a, t, g, c):#psat=0.5, psgc=0.5):
    if sum([a, t, g, c]) != 1:
        raise ArithmeticError('Sum of perentages of A, T, G, and C is not equal to 100 percent!')
    a = ['A' for i in range(round(length * a))]
    t = ['T' for i in range(round(length * t))]
    g = ['G' for i in range(round(length * g))]
    c = ['C' for i in range(round(length * c))]
    unrand = sum([a, t, g, c], [])
    # Help free up memory
    del a, t, g, c
##    at = length * psat
##    gc = length * psgc
##    unrand = sum([['A', 'T'] for i in range(round(at/2))], []) + sum([['G', 'C'] for i in range(round(gc/2))], [])
    seq = []
    for i in range(length):
        rng = (0, len(unrand)-1)
        idx = randint(*rng)
        seq.append(unrand[idx])
        del unrand[idx]
    #print(Seq(''.join(seq)))
    return Seq(''.join(seq))

class gen_proteins(Thread):
    def __init__(self, times, atgc, seqSize):
        Thread.__init__(self)
        self.times = int(times)
        self.atgc = atgc
        self.seqSize = int(seqSize)
        self.start()
    
    def run(self):
        proteins = []
        for i in range(self.times):
            seq = gen_random_proper_seq(self.seqSize, *self.atgc)
            frames = get_frames(seq)
            #for protein in frames:
            #    proteins[len(protein)] = protein
            proteins += [len(protein) for protein in frames]
        self.data = proteins
    pass

def run():
    print('Running...')
    gcrng = (20, 65, 5)
    global threads
    threads = {}
    for gcperc in range(*gcrng):
        at = 100 - gcperc
        at, gc = (at/100, gcperc/100)
        a, t, g, c = (at/2, at/2, gc/2, gc/2)
        threads[gcperc] = gen_proteins(500, [a, t, g, c], 10000)
        print('Thread Started.')
##        #proteins = {}
##        proteins = []
##        for i in range(500):
##            seq = gen_random_proper_seq(10000, a, t, g, c)
##            frames = get_frames(seq)
##            #for protein in frames:
##            #    proteins[len(protein)] = protein
##            proteins += [len(protein) for protein in frames]
##        #max(proteins.keys())
    data = {}
    print('Waiting for threads to end.')
    while threads:
        for key in list(threads.keys()):
            thread = threads[key]
            if not thread.is_alive():
                print('Thread calculating', key, '% GC is Done.')
                data[key] = sum(thread.data) / len(thread.data)
                del threads[key]
    for key in data.keys():
        print('Average Protein Size was', data[key],
              'with', key, '% of random DNA made up of GC.')
    print('Done')

if __name__ == '__main__':
    run()
