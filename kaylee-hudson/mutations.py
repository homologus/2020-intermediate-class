f=open("nucl.fasta")
lines = f.readlines()
seq = ""
extra = len(lines[1]) % 3
print(len(lines[1]))
for i in range(len(lines[1]) - extra):
	seq = seq + lines[1][i]
swap = {'a': 'A', 'g':'G', 't':'T', 'c':'C'}
newSeq = ""
for i in range(len(seq)):
	if(seq[i] == 'a' or seq[i] == 'g' or seq[i] == 'c' or seq[i] == 't'):
		newSeq = newSeq + swap[seq[i]]
	else:
		newSeq = newSeq + seq[i]
proteins = {'TTT' : 'F', 'TTC' : 'F', 'TTA' : 'L', 'TTG' : 'L', 'CTT' : 'L', 'CTC' : 'L', 'CTA' : 'L', 'CTG' : 'L', 'ATT' : 'I', 'ATC' : 'I', 'ATA' : 'I', 'ATG' : 'M', 'GTT' : 'V', 'GTC' : 'V', 'GTA': 'V', 'GTG': 'V', 'TCT' : 'S', 'TCC' : 'S', 'TCA': 'S', 'TCG' : 'S', 'CCT' : 'P', 'CCC' : 'P', 'CCA' : 'P', 'CCG' : 'P', 'ACT' : 'T', 'ACC' : 'T', 'ACA': 'T', 'ACG' : 'T', 'GCT' : 'A', 'GCC' : 'A', 'GCA': 'A', 'GCG' : 'A', 'TAT': 'Y', 'TAC' : 'Y', 'TAA' : '*', 'TAG': '*', 'CAT' : 'H', 'CAC' : 'H', 'CAA' : 'Q', 'CAG' : 'Q', 'AAT' : 'N', 'AAC' : 'N', 'AAA' : 'K', 'AAG' : 'K', 'GAT' : 'D', 'GAC' : 'D', 'GAA' : 'E', 'GAG' : 'E', 'TGT' : 'C', 'TGC' : 'C', 'TGA' : '*', 'TGG' : 'W', 'CGT' : 'R', 'CGC' : 'R', 'CGA' : 'R', 'CGG': 'R', 'AGT' : 'S', 'AGC' : 'S', 'AGA' : 'R', 'AGG' : 'R', 'GGT' : 'G', 'GGC' : 'G', 'GGA' : 'G', 'GGG' : 'G' }
protseq = ""
for i in range(0, len(newSeq), 3):
	protseq = protseq + proteins[newSeq[i:i+3]]
print(protseq)
brandnewseq = ""
for i in range(len(newSeq)):
	if(i==1509):
		brandnewseq = brandnewseq + 'A'
	else:
		brandnewseq = brandnewseq + newSeq[i]
newprotseq = ""
for i in range(0, len(brandnewseq), 3):
	newprotseq = newprotseq + proteins[brandnewseq[i:i+3]]
print(newprotseq)
