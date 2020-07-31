from Bio import SeqIO
reads = list(SeqIO.parse("/share/SARS/covid-patient-CA.fasta", "fasta"))
totalcount = 0
for i in range(0, 3):
	protseq = reads[0].seq[i:].translate()
	aminocount = 0
	newseq = ""
	seqcount = 0
	print("Reading frame " + str(i + 1) + ":")
	for j in range(i, len(protseq)):
		if(protseq[j] == '*'):
			if(aminocount > 100):
				newseq = newseq + protseq[j - aminocount : j + 1]
				seqcount = seqcount + 1
			aminocount = 0
		else:
			aminocount = aminocount + 1
	print(newseq)
	print("Number of protein sequences: " + str(seqcount))
	totalcount = totalcount + seqcount
print("Total number of protein sequences: " + str(totalcount))
