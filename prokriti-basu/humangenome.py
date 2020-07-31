import string
from Bio.Seq import Seq
from Bio import SeqIO

filepath = "/share/Human/chr12.fa"
filepath1 = "/home/prokriti/HumanALDH2.fasta"
filepath2 = "/share/Human/chrX.fa"
filepath3 = "/home/prokriti/greenopsin.fasta"
def find_pos(filepath: str, pos: int):
	reads = list(SeqIO.parse(filepath, "fasta"))
	rf = reads[0].seq
	r = rf[pos]
	return r

def change_pos(filepath: str, pos: int, c: str):
	reads = list(SeqIO.parse(filepath, "fasta"))
	s = ''
	l = len(reads[0])
	rf = list(reads[0].seq)
	rf[pos] =  c
	return s.join(rf)

#a = find_pos(filepath, 111803961)
#print(a)

#b = find_pos(filepath1, 1509)
#print(b)

#c = change_pos(filepath1, 1509,"A")
#print(c)

d = find_pos(filepath2, 154195933)
print(d)

e = find_pos(filepath3, 988)
print(e)

f = change_pos(filepath3, 988, "A")
print(">HumanColorBlind")
print(f)
