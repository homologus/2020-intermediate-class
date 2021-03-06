program = open("/share/SARS/SARS-2020.fasta", "r")
genes = program.readlines()[1:]
genome = []
genome1 = ""

for i in genes:
  i = i.rstrip('\n')
  genome.append(i)
for i in genome:
  for k in i:
    genome1 = genome1+k

codon_table = {"TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L", "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M", "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V", "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P", "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "GCG": "A", "GCT": "A", "GCC": "A", "GCA": "A", "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*", "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q", "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E", "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W", "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R", "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R", "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G"}
protein_chain1 = ""
protein_chain2 = ""
protein_chain3 = ""
reading_frame1 = []
reading_frame2 = []
reading_frame3 = []

for i in range(0, len(genome1), 3):
  reading_frame1.append(genome1[i:i+3])
for i in range(1, len(genome1), 3):
  reading_frame2.append(genome1[i:i+3])
for i in range(2, len(genome1), 3):
  reading_frame3.append(genome1[i:i+3])

for i in reading_frame1:
  if i in codon_table:
    protein_chain1 = protein_chain1 + codon_table[i]
for i in reading_frame2:
  if i in codon_table:
    protein_chain2 = protein_chain2 + codon_table[i]
for i in reading_frame3:
  if i in codon_table:
    protein_chain3 = protein_chain3 + codon_table[i]

answer = 0
long_protein_sequences = []
protein_chain1 = protein_chain1.split('*')
protein_chain2 = protein_chain2.split('*')
protein_chain3 = protein_chain3.split('*')
seq_length = 0
begin = []
end = []

for i in protein_chain1:
  seq_length = seq_length+len(i)
  if len(i)>99:
    answer = answer+1
    long_protein_sequences.append(i)
    begin.append(seq_length)
    end_of_seq = seq_length+len(i)+1
    end_of_seq = end_of_seq * 3
    end.append(end_of_seq)
seq_length = 0
for i in protein_chain2:
  seq_length = seq_length+len(i)
  if len(i)>99:
    answer = answer+1
    long_protein_sequences.append(i)
    begin.append(seq_length)
    end_of_seq = seq_length+len(i)+1
    end_of_seq = 1+end_of_seq * 3
    end.append(end_of_seq)
seq_length = 0
for i in protein_chain3:
  seq_length = seq_length+len(i)
  if len(i)>99:
    answer = answer+1
    long_protein_sequences.append(i)
    begin.append(seq_length)
    end_of_seq = seq_length+len(i)+1
    end_of_seq = 2+end_of_seq * 3
    end.append(end_of_seq)

print ("There are " + str(answer) + " long protein sequences.")
print ("The protein sequences are: ", long_protein_sequences)

reverse_codon_table = {"F": "TTT", "F": "TTC", "L": "TTA", "L": "TTG", "L": "CTT", "L": "CTC", "L": "CTA", "L": "CTG", "I": "ATT", "I": "ATC", "I": "ATA", "M": "ATG", "V": "GTT", "V": "GTC", "V": "GTA", "V": "GTG", "S": "TCT", "S": "TCC", "S": "TCA", "S": "TCG", "P": "CCT", "P": "CCC", "P": "CCA", "P": "CCG", "T": "ACT", "T": "ACC", "T": "ACA", "T": "ACG", "A": "GCG", "A": "GCT", "A": "GCC", "A": "GCA", "Y": "TAT", "Y": "TAC", "*": "TAA", "*": "TAG", "H": "CAT", "H": "CAC", "Q": "CAA", "Q": "CAG", "N": "AAT", "N": "AAC", "K": "AAA", "K": "AAG", "D": "GAT", "D": "GAC", "E": "GAA", "E": "GAG", "C": "TGT", "C": "TGC", "*": "TGA", "W": "TGG", "R": "CGT", "R": "CGC", "R": "CGA", "R": "CGG", "S": "AGT", "S": "AGC", "R": "AGA", "R": "AGG", "G": "GGT", "G": "GGC", "G": "GGA", "G": "GGG"}

print (begin)
print (end)
