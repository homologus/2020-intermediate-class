from Bio.Seq import Seq
from Bio import SeqIO

#finding the nucleotide at position 111803962, which determines the asian flush.
reads = list(SeqIO.parse("/share/Human/chr12.fa", "fasta"))

print(reads[0].seq[111803961:111803962])

#modifying the ALDH2 gene to make it have asian flush
#creating file and storing new ALDH2 gene
file = open("nuclu.fasta", 'w+')
gene = list(SeqIO.parse("nucl.fasta", "fasta"))

#overwriting "G" --> "A" at position 1509
new = str(gene[0].seq[:1509]) +"A"+ str(gene[0].seq[1510:])
#adding sequence to new fasta file "nuclu.fasta"
file.write(">Mutated ALDH2 gene \n")
file.write(new)
file.close()
