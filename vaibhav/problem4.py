program = open("/share/Ecoli/GCA_000005845.2_ASM584v2_genomic.fna", "r")
genes = program.readlines()[1:]
genome = []
genome1 = []
for i in genes:
  i = i.rstrip('\n')
  genome.append(i)
for i in genome:
  genome1.append(i)
genome2 = ""

for i in range(len(genome1)):
  if i == "G":
    if i+1 == "A":
      if i+2 == "C":
        if i+8 == "G":
          if i+9 == "T":
            if i+10 == "C":
              genome2 = genome2 + genome1[i:i+6] + " " + genome1[i+6, i+11]
              print (genome2)
  else:
    genome2 = genome2+genome1[i]

genome3 = genome2.split(" ")
lengths = []

for i in genome3:
  lengths.append(len(i))

#print (genome3)
print (lengths)
