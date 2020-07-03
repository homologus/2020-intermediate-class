from Bio.Seq import Seq
from Bio import SeqIO

def getNum( e ):
    return int(e[0] )

with open("matches.txt", "r") as f:
    match = f.read()

match = match[match.index(">") : match.index("Lambda")]

test = match.split("Score")

test.pop(0)

table = [] 
for key in test:
    key = key[key.index("Query") :].strip()
    exon = key.split("\n")
    
    l1 = exon[0].split(" ")
    l2 = exon[2].split(" ")
    l3 = exon[len(exon) - 3].split(" ")
    l4 = exon[len(exon) - 1].split(" ")
    
    table.append( (int(l1[2]),  int(l3[len(l3) - 1]),int ( l2[2]),  int(l4[len(l4) - 1]) ) )


table.sort(key = getNum)



search =  111803962 #111803962

searchLoc = 0
for i, k in enumerate(table):
    if k[3] >= search and k[2] <= search:
        searchLoc = i
        break

search = table[searchLoc][1] - (table[searchLoc][3] - search) - 1

print(search)




ALDH2 = list(SeqIO.parse("nucl.fasta", "fasta"))
ALDH2 = ALDH2[0].seq.lower()
proNoMutALDH2 = ALDH2.translate()



print(ALDH2[search])
ALDH2 = ALDH2[0 : search] + "a" + ALDH2[search + 1:]

proYesMutALDH2 = ALDH2.translate()
print(proNoMutALDH2)
print(proYesMutALDH2)

for i in range(0, len(proNoMutALDH2), 1):
    if proNoMutALDH2[i : i + 1] == proYesMutALDH2[i : i + 1]:
         count = 0;
    else:
        print(proNoMutALDH2[i : i + 1], "    ", proYesMutALDH2[i : i + 1]) 

