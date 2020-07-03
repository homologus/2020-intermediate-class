
from Bio.Seq import Seq


from Bio import SeqIO


fastaFile = list(SeqIO.parse("/share/Human/chr12.fa", "fasta"))
seq = fastaFile[0].seq

print(seq[111803962 - 1])


