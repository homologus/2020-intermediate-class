from Bio import SeqIO 

print('SARS-2020')
covid = list(SeqIO.parse("/share/SARS/SARS-2020.fasta","fasta"))
l = []

sequence = covid[0].seq
s = sequence.translate()
frameOne = s.rsplit('*')
sequence = covid[0].seq[1:]
ss = sequence.translate()
frameTwo = ss.rsplit('*')
sequence = covid[0].seq[2:]
sss = sequence.translate()
frameThree = sss.rsplit('*')
for i in range(len(frameOne)):
	if len(frameOne[i]) > 100:		
		l.append(frameOne[i])
		hi = str(frameOne[i])
		print(s.index(hi))
		print('-')
for i in range(len(frameTwo)):
	if len(frameTwo[i]) > 100:
		l.append(frameTwo[i])		
		print(frameTwo[i])
		print('-')
for i in range(len(frameThree)):
	if len(frameThree[i]) > 100:
		l.append(frameThree[i])
		print(frameThree[i])
		print('-')
print(" ")
print('SARS Bat')
covid = list(SeqIO.parse("/share/SARS/hku-bat.fasta","fasta"))

sequence = covid[0].seq
frameOne = sequence.translate()
frameOne = frameOne.rsplit('*')
sequence = covid[0].seq[1:]
frameTwo = sequence.translate()
frameTwo = frameTwo.rsplit('*')
sequence = covid[0].seq[2:]
frameThree = sequence.translate()
frameThree = frameThree.rsplit('*')
for i in range(len(frameOne)):
	if len(frameOne[i]) > 100:
		l.append(frameOne[i])
		print(frameOne[i])
		print('-')
for i in range(len(frameTwo)):
	if len(frameTwo[i]) > 100:
		l.append(frameTwo[i])
		print(frameTwo[i])
		print('-')
for i in range(len(frameThree)):
	if len(frameThree[i]) > 100:
		l.append(frameThree[i])
		print(frameThree[i])
		print('-')
print(" ")
print('SARS Pangolin')
covid = list(SeqIO.parse("/share/SARS/pangolin-cov.fasta","fasta"))

sequence = covid[0].seq
frameOne = sequence.translate()
frameOne = frameOne.rsplit('*')
sequence = covid[0].seq[1:]
frameTwo = sequence.translate()
frameTwo = frameTwo.rsplit('*')
sequence = covid[0].seq[2:]
frameThree = sequence.translate()
frameThree = frameThree.rsplit('*')
for i in range(len(frameOne)):
	if len(frameOne[i]) > 100:
		l.append(frameOne[i])
		print(frameOne[i])
		print('-')
for i in range(len(frameTwo)):
	if len(frameTwo[i]) > 100:
		l.append(frameTwo[i])
		print(frameTwo[i])
		print('-')
for i in range(len(frameThree)):
	if len(frameThree[i]) > 100:
		l.append(frameThree[i])
		print(frameThree[i])
		print('-')
print(" ")
print('SARS MERS')
covid = list(SeqIO.parse("/share/SARS/mers.fasta","fasta"))

sequence = covid[0].seq
frameOne = sequence.translate()
frameOne = frameOne.rsplit('*')
sequence = covid[0].seq[1:]
frameTwo = sequence.translate()
frameTwo = frameTwo.rsplit('*')
sequence = covid[0].seq[2:]
frameThree = sequence.translate()
frameThree = frameThree.rsplit('*')
for i in range(len(frameOne)):
	if len(frameOne[i]) > 100:
		l.append(frameOne[i])
		print(frameOne[i])
		print('-')
for i in range(len(frameTwo)):
	if len(frameTwo[i]) > 100:
		l.append(frameTwo[i])
		print(frameTwo[i])
		print('-')
for i in range(len(frameThree)):
	if len(frameThree[i]) > 100:
		l.append(frameThree[i])
		print(frameThree[i])
		print('-')


