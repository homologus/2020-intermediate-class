program = open('/share/SARS/SARS-2020.fasta', 'r')
genes = program.readlines()[1:]
genome = []
for i in genes:
  i = i.rstrip('\n')
  genome.append(i)

a=0
c=0
g=0
t=0

for i in genome:
  for j in i:

    if j == "A":
      a = a+1

    if j == "C":
      c = c+1

    if j == "G":
      g = g+1

    if j == "T":
      t = t+1

print ("A: "+str(a))
print ("C: "+str(c))
print ("G: "+str(g))
print ("T: "+str(t))

gc = g + c
total = a+c+g+t
acgt = gc/total *100

print ("The GC Content is " + str(acgt) + "%.")
