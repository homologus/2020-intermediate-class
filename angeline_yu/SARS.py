from Bio import SeqIO
from Bio.Seq import Seq

#Accessing SARS Genome and storing it in reads
reads = list(SeqIO.parse("/share/SARS/SARS-2020.fasta", "fasta"))

#translating the 3 reading frames
translation1 = reads[0].seq.translate()
translation2 = reads[0].seq.translate()
translation3 = reads[0].seq.translate()
for i in reads:
        #reading frame 1
        translation1 = i.seq[0:].translate()
        #reading frame 2
        translation2 = i.seq[1:].translate()
        #reading frame 3
        translation3 = i.seq[2:].translate()

#total genes
total = 0
#count variable acts as place holder to check if sequence is 100 long
count = 0
#toReturn variable keeps track of how many 100 amino acid sequences there are
toReturn = 0
#saves previous i location
prev_i = 0
#For loops checks for end codens ("*")
#reading frame 1
for i in range(len(translation1)):
        if(str(translation1[i]) != "*"):
                count += 1
                continue
        else:
                if(count >= 100):
                        toReturn += 1
                        print("\n>Amino Acid Sequence " +str(toReturn))
                        print("Start Location = " +str(prev_i))
                        print("End Location = " +str(i+1))
                        print("Length = " +str((i+1)-prev_i))
                        print(translation1[prev_i+1:i]+"\n")
                prev_i = i
                count = 0

print("\nThere are " + str(toReturn) + " 100 amino acid sequences in the SARS-CoV2 Genome in reading frame 1. \n")
total += toReturn

count = 0
toReturn = 0
prev_i = 0
for i in range(len(translation2)):
        if(str(translation2[i]) != "*"):
                count += 1
                continue
        else:
                if(count >= 100):
                        toReturn += 1
                        print("\n>Amino Acid Sequence " +str(toReturn))
                        print("Start Location = " +str(prev_i))
                        print("End Location = " +str(i+1))
                        print("Length = " +str((i+1)-prev_i))
                        print(translation2[prev_i+1:i]+"\n")
                prev_i = i
                count = 0

print("There are " + str(toReturn) + " 100 amino acid sequences in the SARS-CoV2 Genome in reading frame 2. \n")
total += toReturn

count = 0
toReturn = 0
prev_i = 0
for i in range(len(translation3)):
        if(str(translation3[i]) != "*"):
                count += 1
                continue
        else:
                if(count >= 100):
                        toReturn += 1
                        print("\n>Amino Acid Sequence " +str(toReturn))
                        print("Start Location = " +str(prev_i))
                        print("End Location = " +str(i+1))
                        print("Length = " +str((i+1)-prev_i))
                        print(translation3[prev_i+1:i]+"\n")
                prev_i = i
                count = 0

print("There are " + str(toReturn) + " 100 amino acid sequences in the SARS-CoV2 Genome in reading frame 3. \n")
total += toReturn

print("There are " +str(total) + " total 100 amino acid sequences in the SARS-CoV2 Genome. \n")

