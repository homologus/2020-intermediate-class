seq=open("seq.fasta")
lines=seq.readlines()
length=len(lines)-1
length1=len(lines[1])-2
lastlinelen=len(lines[-1])
print(length*length1+lastlinelen)
