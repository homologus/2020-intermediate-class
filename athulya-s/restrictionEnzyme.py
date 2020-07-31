import re

file = open('/share/Ecoli/GCA_000005845.2_ASM584v2_genomic.fna','r')
seq = file.readlines()

sequence = ''
count = 0
for i in range(1,len(seq)):
	sequence = sequence + seq[i].rstrip('\n')

a = re.split('GAC[ACGT][ACGT][ACGT][ACGT][ACGT]GTC',sequence)
print('Number of fragments split: ' + str(len(a)-1))
for i in range(len(a)):
	print('Length of each fragment:' + str(len(a[i])))

