import re
fileName = open('/share/SARS/SARS-2020.fasta','r')
seq = fileName.readlines()

nucleotides = {'A':0,'T':0,'C':0,'G':0}
count = 0

for i in range(len(seq)):
        count = len(seq[i])
        for j in range(count):
                if nucleotides.get(seq[i][j])!= None:
                        nucleotides[seq[i][j]] = nucleotides[seq[i][j]]+1
print(nucleotides)
print(((nucleotides['C']+nucleotides['G'])/(nucleotides['G']+nucleotides['C']+nucleotides['A']+nucleotides['T']))*100),print('%')
print(len(seq))
