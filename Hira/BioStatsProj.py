import random
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqUtils import GC

def generateSequence(gc):
    g = gc*50
    c = g
    a = (100 - g*2)/2
    t = a
    seq = ""
    for i in range(0, 10000, 1):
        num = random.randint(1, 100)
        if num <= g:
            seq = seq + 'g'
        elif num <= g + c:
            seq = seq + 'c'       
        elif num <= g + c + a:
            seq = seq + 'a'
        else:
            seq = seq + 't'
    return Seq(seq)


def findLongestProtein(seq):
    i = 0
    proteins = [""]
    while i <= len(seq)- len(seq) % 3:
        protein = seq[i : ].translate(to_stop = True)
        if len(protein) > len(proteins[ len(proteins) - 1] ):
            proteins.append(protein)
            i += len(protein) * 3
        else:
            i += 1
    return proteins[ len(proteins) - 1]


for i in range(15, 66, 5):
    seq = generateSequence(i/100)
    protein = findLongestProtein(seq)
    print("GC: ~", i, "% (", GC(seq) ,"%),  Length: ", len(protein))


fastaFile  = list(SeqIO.parse("/share/Plasmodium/PfGB4.Jul2015.fasta", "fasta"))
seq1 = fastaFile[0].seq.lower()



fastaFile = list(SeqIO.parse("/share/Ecoli/GCA_000005845.2_ASM584v2_genomic.fna", "fasta"))
seq2 = fastaFile[0].seq.lower()


print("Plasmodium GC%: ", GC(seq1), "%")

print("Ecoli GC% : ", GC(seq2), "%")


