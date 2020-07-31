import random
from Bio import SeqIO
from Bio.Seq import Seq

reads = list(SeqIO.parse("/share/SARS/hku-bat.fasta", "fasta"))

dna_sequence1 = reads[0].seq[0:]
dna_sequence2 = reads[0].seq[1:]
dna_sequence3 = reads[0].seq[2:]
protein1 = dna_sequence1.translate()
protein2 = dna_sequence2.translate()
protein3 = dna_sequence3.translate()

protein1 = protein1.split("*")
protein2 = protein2.split("*")
protein3 = protein3.split("*")

answer = 0
protein_chain = []
for i in protein1:
  if len(i) > 99:
    print(i)
    print("")
    answer = answer + 1

for i in protein2:
  if len(i) > 99:
    print(i)
    print("")
    answer = answer + 1

for i in protein3:
  if len(i) > 99:
    print(i)
    print("")
    answer = answer + 1

print("There are " + str(answer) + " protein sequences.")
#print(protein_chain)

