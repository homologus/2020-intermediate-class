from Bio import SeqIO
from Bio.Seq import Seq

#Accessing SARS Genome and storing it in reads
reads = list(SeqIO.parse("/share/SARS/SARS-2020.fasta", "fasta"))
read = reads[0]
print(read.seq.reverse_complement())
print(reads[0].seq.reverse_complement())

#translating the 3 reading frames
translation1 = reads[0].seq.translate()
translation2 = reads[0].seq.translate()
translation3 = reads[0].seq.translate()

#reverse translating the 3 reading frames
rTranslation1 = reads[0].seq.reverse_complement().translate()
rTranslation2 = reads[0].seq.reverse_complement().translate()
rTranslation3 = reads[0].seq.reverse_complement().translate()

for i in reads:
	#reading frame 1
	translation1 = i.seq[0:].translate()
	rTranslation1 = i.seq[0:].reverse_complement().translate()
	#reading frame 2
	translation2 = i.seq[1:].translate()
	rTranslation2 = i.seq[1:].reverse_complement().translate()
	#reading frame 3
	translation3 = i.seq[2:].translate()
	rTranslation3 = i.seq[2:].reverse_complement().translate()

sequences = [translation1, translation2, translation3]
rSequences = [rTranslation1, rTranslation2, rTranslation3]

#print("1: \n"+rTranslation1)
#print("2: \n"+rTranslation2)
#print("3: \n"+rTranslation3)

print("\nRegular Translation Result: \n")
#Regular Translation
#total genes
total = 0
for sequence in sequences:
	#count variable acts as place holder to check if sequence is 100 long
	count = 0
	#toReturn variable keeps track of how many 100 amino acid sequences there are
	toReturn = 0
	#saves previous i location
	prev_i = 0
	#reading frames
	for i in range(len(sequence)):
		#For loops checks for end codens ("*")
		if(str(sequence[i]) != "*"):
			count += 1
			continue
		#If end coden, stop and check if length over 100
		else:
			if(count >= 100):
				toReturn += 1
				print("\n>Amino Acid Sequence " +str(toReturn))
				print("Start Location = " +str(prev_i))
				print("End Location = " +str(i+1))
				print("Length = " +str((i+1)-prev_i))
				print(sequence[prev_i+1:i]+"\n")
			prev_i = i
			count = 0

	print("\nThere are " + str(toReturn) + " 100 amino acid sequences in the SARS-CoV2 Genome " 
		+ "in reading frame " +str(sequences.index(sequence)+1) +". \n")
	total += toReturn


print("There are " +str(total) + " total 100 amino acid sequences in the SARS-CoV2 Genome. \n")

print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print("Reverse Translation Result: \n")
#Reverse Translation
#total genes
total = 0
for sequence in rSequences: 
	#count variable acts as place holder to check if sequence is 100 long
	count = 0
	#toReturn variable keeps track of how many 100 amino acid sequences there are
	toReturn = 0
	#saves previous i location
	prev_i = 0
	#reading frames
	for i in range(len(sequence)):
		#For loops checks for end codens ("*")
		if(str(sequence[i]) != "*"):
			count += 1
			continue
		#If end coden, stop and check if length over 100
		else:
			if(count >= 100):
				toReturn += 1
				print("\n>Amino Acid Sequence " +str(toReturn))
				print("Start Location = " +str(prev_i))
				print("End Location = " +str(i+1))
				print("Length = " +str((i+1)-prev_i))
				print(sequence[prev_i+1:i]+"\n")
			prev_i = i
			count = 0

	print("\nThere are " + str(toReturn) + " 100 amino acid sequences in the SARS-CoV2 Genome "
		+ "in reading frame " +str(rSequences.index(sequence)+1) +". \n")
	total += toReturn

print("There are " +str(total) + " total 100 amino acid sequences in the SARS-CoV2 Genome. \n")

