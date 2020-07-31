from Bio import SeqIO
from Bio.Seq import Seq
reads = list(SeqIO.parse("/share/Human/chr12.fa", "fasta"))
read = list(SeqIO.parse('nucl.fasta','fasta'))

s = ''
for i in range(len(read[0].seq)):
	s = s + read[0][i]
a = s[:1509] + 'A' + s[1510:]

readds = Seq(a)
sa = Seq(s)
org = sa.translate()
mutation = readds.translate()

print(org)
print(mutation)
for i in range(len(org)):
	if org[i] != mutation[i]:
		print('mutation')
		print(org[i])
		print(mutation[i])
		print('-')

