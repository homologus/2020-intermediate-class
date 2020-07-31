from Bio import SeqIO
from Bio.Seq import Seq

#This function will take one nucleotide sequence and translate it
def amino_acids(virus):
	#initiazing the 3 reading frames
	translation1 = virus[0].seq.translate()
	translation2 = virus[0].seq.translate()
	translation3 = virus[0].seq.translate()

	for i in virus:
		#reading frame 1
		translation1 = i.seq[0:].translate()
		#reading frame 2
		translation2 = i.seq[1:].translate()
		#reading frame 3
		translation3 = i.seq[2:].translate()

	sequences = [translation1, translation2, translation3]

	print("\nRegular Translation Result: ")
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
					print("Reading Frame: "+str(sequences.index(sequence)+1))
					print("Start Location = " +str(prev_i))
					print("End Location = " +str(i+1))
					print("Length = " +str((i+1)-prev_i))
					print(sequence[prev_i+1:i]+"\n")

				prev_i = i
				count = 0

		print("\nThere are " + str(toReturn) + " 100 amino acid sequences in the SARS-CoV2 Genome "
			+ "in reading frame " +str(sequences.index(sequence)+1) +". \n")
		total += toReturn


	print("There are " +str(total) + " total 100 amino acid sequences in the SARS-CoV2 Genome. ")


#List of different Sequences
#This list will hold all the original sequences
reads = []
#Beta Coronavirus
reads.append(list(SeqIO.parse("/share/SARS/betacov.fasta", "fasta")))
#SARS Coronavirus from 2003
reads.append(list(SeqIO.parse("/share/SARS/SARS-urbani.fasta", "fasta")))
#MERS Coronavirus
reads.append(list(SeqIO.parse("/share/SARS/mers.fasta", "fasta")))
#Bat Coronoavirus
reads.append(list(SeqIO.parse("/share/SARS/hku-bat.fasta", "fasta")))
#Bat Coronavirus 2
reads.append(list(SeqIO.parse("/share/SARS/nature-2013.fasta", "fasta")))
#Coronavirus from pangolin
reads.append(list(SeqIO.parse("/share/SARS/pangolin-cov.fasta", "fasta")))
#SARS-2 Coronavirus from Wuhan
reads.append(list(SeqIO.parse("/share/SARS/SARS-2020.fasta", "fasta")))
#SARS-2 Coronavirus from first patient in CA
reads.append(list(SeqIO.parse("/share/SARS/covid-patient-CA.fasta", "fasta")))

#Coronavirus names to be displayed
names = [
	"Beta Coronavirus", "SARS Coronavirus from 2003", "MERS Coronavirus",
        "Bat Coronavirus", "Bat Coronavirus 2", "Coronavirus from pangolin",
        "SARS-2 Coronavirus from Wuhan", "SARS-2 Coronavirus from first patient in CA"
	]
#This Variable will keep track of which name to display
name = 0
#For loop to repeat for every virus sequences
for virus in reads:
	print("\n***********************************************************************************")
	print(names[name]+ ": ")
	name += 1
	amino_acids(virus)

print("***********************************************************************************")
