#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Gets specific nucleotides in DNA from a file.
# Programmed by Samuel Davenport

from Bio.Seq import Seq
from Bio import SeqIO
from random import randint, choice
import os
#from threading import Thread

NAME = 'Recoupler'
__version__ = '0.0.1'

#DATA = [i[1:] if not ' ' in i[1:] else i[1:][:i.index(' ')] for i in DATA.split(',')]

FILE = '/share/Human/chr12.fa'
FTYPE = 'fasta'

progs = ['blastn', 'blastp', 'blastx', 'tblastn', 'tblastx']

FILE = 'SARS-2020.fasta'

if os.sys.argv[1:]:
    if '-f' in os.sys.argv:
        idx = os.sys.argv.index('-f')
        FILE = os.sys.argv[idx+1]
    if '-t' in os.sys.argv:
        idx = os.sys.argv.index('-t')
        FTYPE = os.sys.argv[idx+1]

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

def recouple(lst, scanLen=8, start=''):
    """Program that will re-couple text from a list that has been broken into many overlapping parts."""
    # Get all the seperate chuncks of data in a list, no repeats.
    data = []
    [data.append(i) for i in lst if not i in data]
    if scanLen < 4:
        scanLen = min(round(min([len(i) for i in data])/2), 4)
        print('scanLen =',scanLen)
    # Get all the words into a dictionary
    wordsLst = {}
    beginsings = {}
    # For each thing to test in the dictionary
    for test in data:
        # If we haven't seen this word before, set it as a possible
        # beginning value
        if not test[0:scanLen] in beginsings.keys():
            beginsings[test[0:scanLen]] = True
        # For each index point our test string,
        for i in range(len(test)-scanLen):
            # Get the word it could be
            word = test[i:i+scanLen]
            # Add the word to the words dict, pointing to the last letter of
            # our test string.
            wordsLst[word] = test[i+scanLen]
            # If it's not the beginning
            if i > 0:
                beginsings[word] = False
    # If the starting word was not defined,
    if not start:
        # Set it to a random valid word
        for word in beginsings.keys():
            if beginsings[word]:
                test = str(word)
                break
    else:
        # Otherwise, try to use it
        test = str(start)
        cstart = False
        yes = False
        # If it's invalid, break
        if not test in wordsLst.keys():
            for word in wordsLst.keys():
                if test in word:
                    test = word
                    yes = True
                    break
        if not yes:
            raise ValueError('Start word argument is invalid, does not exist in word list with scanLen. Try having start word len be equal to scanLen argument number.')
    # Free up some memory
    del beginsings
    # Get data that points to eachother
    string = test
    count = len(wordsLst.keys())*scanLen
    while count:
        # If the test word is valid
        if test in wordsLst.keys():
            # Add it to the string
            string += wordsLst[test]
            test = test[1:scanLen]+wordsLst[test]
            count -= 1
        else:# Otherwise, decrement scan by one
            count = 0
##    # Once we have data that points to other bits, put them together properly
##    data = str(string[0])
##    # For each bit of data to add,
##    for i in string[1:]:
##        # For all possible index positions in the new bit of data,
##        for ii in range(len(i)):
##            # Get the current data's tail with ii
##            piece = data[-ii:]
##            # If the tail of this index is the start of this bit of data,
##            if piece == i[:len(piece)]:
##                # Combine the two chucks properly and stop checking indexes for this new chunk.
##                data = data[:len(data)-ii] + i
##                break
    return string

def getNucleotide(seq, position):
    if position < len(seq):
        return seq[position-1]
    return None

def run():
    print('Running...')
    dna_seq = SeqIO.read(FILE, FTYPE)
    dna_seq = dna_seq.seq[:-(len(dna_seq) % 3)]
    while True:
        position = input('position : ')
        if not position.isnumeric():
            break
        change = input('change : ')
        one = dna_seq.translate()
        seq2 = str(dna_seq)
        seq2[position] = change
        two = Seq(seq2).translate()
        change -= change % 3
        idx = int(change / 3)
        print(one[idx], two[idx])
        #print(getNucleotide(dna_seq, int(position)))
    print('Done')

if __name__ == '__main__':
    run()
