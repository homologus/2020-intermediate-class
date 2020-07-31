from Bio import SeqIO
reads = list(SeqIO.parse("/share/Human/chr12.fa", "fasta"))
print(reads[0].seq[111803961])
