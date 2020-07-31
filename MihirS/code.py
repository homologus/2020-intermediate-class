from Bio import SeqIO
from Bio.Seq import Seq
SARS = list(SeqIO.parse("small-genome", "fasta"))

count = 0
SARS = SARS[0]
print("reading frame 1")
protein = SARS.translate()
print (protein)
for i in protein:
	if i != '*':
		count = count + 1
	elif i == '*':
		if count >= 100:
			print (protein[i-count:i])
		count = 0

count = 0
SARS = SARS[1:len(SARS)-2]
print ("reading frame 2")
protein = SARS.translate()
print (protein)
for i in protein:
	if i != '*':
		count = count + 1
	elif i == '*':
		if count >= 100:
			print (protein[i-count:i])
		count = 0

count = 0
SARS = SARS[1:len(SARS)-1]
print ("reading frame 3")
print(SARS.translate())
