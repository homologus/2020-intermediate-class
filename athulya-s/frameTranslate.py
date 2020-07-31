file = open('/share/SARS/SARS-2020.fasta','r')
seq = file.readlines()
proteins = {'TTT':'F','TTC':'F','TTA':'L','TTG':'L','TCT':'S','TCC':'S','TCA':'S','TCG':'S','TAT':'Y','TAC':'Y','GTC':'V','TAA':'*','TAG':'*','TGT':'C','TGC':'C','TGA':'*','TGG':'W','CTT':'L','CTC':'L','CTA':'L','CTG':'L','CCT':'P','CCC':'P','CCA':'P','CCG':'P','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q','CGT':'R','CGC':'R','CGA':'R','CGG':'R','ATT':'I','ATC':'I','ATA':'I','ATG':'M','ACT':'T','ACC':'T','ACA':'T','ACG':'T','AAT':'N','AAC':'N','AAA':'K','AAG':'K','AGT':'S','AGC':'S','AGA':'R','AGG':'R','GTT':'V','GTG':'V','GTA':'V','GTG':'V','GCT':'A','GCC':'A','GCA':'A','GCG':'A','GAT':'D','GAC':'D','GAA':'E','GAG':'E','GGT':'G','GGC':'G','GGC':'G','GGA':'G','GGG':'G'}

sequence = ""
pro = ''

#turn seq into a string -> sequence
for i in range(1,len(seq)):
	length = len(seq[i])
	seq[i] = seq[i].rstrip('\n')
	sequence = sequence + seq[i]

#frame translate @ 0
for j in range(0,len(sequence), 3):
	if len(sequence[j:j+3]) == 3:
		pro = pro + proteins[sequence[j:j+3]]
p = pro.split('*')
for k in range(len(p)):
	if len(p[k]) > 100:
		print(str(pro.index(p[k])) +' '+ str(pro.index(p[k])+len(p[k])))	
pro = ''
print("next frame")

#frame translate @ 1
for j in range(1, len(sequence), 3):
	if len(sequence[j:j+3]) == 3:
		pro = pro + proteins[sequence[j:j+3]]
p = pro.split('*')
for k in range(len(p)):
	if len(p[k]) > 100:
		print( str ( pro.index ( p[k] ) )+ ' ' + str ( pro.index ( p[k] )+len ( p[k] ) ) )
pro = ''
print("next frame")

#frame translate @ 2
for j in range(2, len(sequence), 3):
	if len(sequence[j:j+3])==3:
		pro = pro + proteins[sequence[j:j+3]]
p = pro.split('*')
for k in range(len(p)):
	if len(p[k]) > 100:
		print(str(pro.index(p[k])) + ' ' + str(pro.index(p[k])+len(p[k])))
