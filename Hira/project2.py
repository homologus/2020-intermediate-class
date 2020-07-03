
from Bio.Seq import Seq


from Bio import SeqIO

dna = Seq("GATCGTAGATAGTGCGCGCGTAGAGGAGAGATAGAGAGAGGAGATAGAGAT")


fastaFile = list(SeqIO.parse("/share/SARS/seq.fasta", "fasta"))
seq = fastaFile[0].seq

proteins = []



i = 0
while i <= len(seq) - 300:
    protein = seq[i : ].translate(to_stop = True) 
    if len(protein) >= 100:
        proteins.append(protein)
        i += len(protein) * 3
        print(protein+ "\n")
    else:
        i += 1;


#for i in range(0, len(seq) -300, 1):
#    protein = seq[i : ].translate(to_stop = True)
#    if len(protein) >= 100:
#        proteins.append(protein)
#        print(protein+ "\n")
