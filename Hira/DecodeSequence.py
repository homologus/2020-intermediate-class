
from Bio.Seq import Seq


from Bio import SeqIO



def findNextStr(fragments, txt):
    for i, key in enumerate(fragments):
        if key[0 : len(txt) ] == txt:
            return key
    return None

         
def findBeforeStr(fragment, txt):
    for i, key in enumerate(fragments):
        key = key.strip()
        if key[len(key) -  len(txt) : ] == txt:
            print(txt) 
            return key.strip()
    return None


with open("Sequence.txt", "r") as f:
    fragments =  f.read()

fragments =  fragments.replace("\n", ",")
fragments = fragments.split(",")

sequence = fragments[0]; 
x = 0; #the start value


while len(fragments) > 0:
    partStr = fragments[x]
    i = 3
    nextStr = None
    fragments.pop(x)
    while nextStr == None and i < len(partStr) :
        nextStr = findNextStr( fragments, partStr[len(partStr) - i : ])
        i += 1
    if i >= len(partStr):
        break; 
    x = fragments.index(nextStr)
    sequence = sequence + nextStr[i - 1: ]
#    print(sequence)
    
sequence = sequence.strip()
while len(fragments) > 0:
    i = 2
    beforeStr = None
#    fragments.pop(x)
    while beforeStr == None and i < 15 :
#        beforeStr = findBeforeStr( fragments, partStr[0 : i ])
        beforeStr = findBeforeStr( fragments, sequence[0 : i ])
        print(i)
        i += 1
    if i >= 15:
        print("herererere?????")
        break;
    x = fragments.index(beforeStr)
    fragments.pop(x)
    sequence =  beforeStr.strip() +  sequence[i - 1 :]
    print("|", sequence, "  ", i - 1, " ")
  
