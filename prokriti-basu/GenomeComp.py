import string
from Bio.Seq import Seq
from Bio import SeqIO

filepath = "/home/prokriti/HumanALDH2.fasta"
filepath1 = "/home/prokriti/MutatedHumanALDH2.fasta"
filepath2 = "/home/prokriti/mutatedgreenopsin.fasta"
filepath3 = "/home/prokriti/greenopsin.fasta"
def read(filepath: str):
        reads = list(SeqIO.parse(filepath, "fasta"))
        for i in reads:
                rf1 = i.seq[0:]
        return rf1

def  protein(s):
        p = s.transcribe()
        p = p.translate()
        return p

#a = read(filepath)
#b = read(filepath1)
#c = protein(a)
#d = protein(b)
#print("Non Asian Flush")
#print(c)
#print("Asian Flush")
#print(d)

r1 = read(filepath2)
r2 = read(filepath3)
p1 = protein(r1)
p2 = protein(r2)
print(">HumanGreenOpsin")
print(p2)
print(">MutatedGreenOpsin")
print(p1)
