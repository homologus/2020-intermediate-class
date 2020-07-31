from Bio import SeqIO
reads = list(SeqIO.parse("/share/SARS/SARS-2020.fasta", "fasta"))
print (reads)

reads1 = reads[0].seq[0:]
reads2 = reads[0].seq[1:]
reads3 = reads[0].seq[2:]
reads1 = reads1.translate()
reads2 = reads2.translate()
reads3 = reads3.translate()

reads1 = reads1.split("*")
reads2 = reads2.split("*")
reads3 = reads3.split("*")
answer = 0
long_protein_sequences = []

for i in reads1:
  if len(i)>99:
    answer = answer+1
    long_protein_sequences.append(i)

for i in reads2:
  if len(i)>99:
    answer = answer+1
    long_protein_sequences.append(i)

for i in reads3:
  if len(i)>99:
    answer = answer+1
    long_protein_sequences.append(i)

print ("There are " + str(answer) + " long protein sequences.")
print ("The protein sequences are: ", long_protein_sequences)

