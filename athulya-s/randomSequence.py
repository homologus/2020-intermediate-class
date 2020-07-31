from Bio import SeqIO
import random
from Bio.SeqUtils import GC 
from Bio.Seq import Seq
covid = list(SeqIO.parse("/share/SARS/SARS-2020.fasta","fasta"))

gc = GC(covid[0].seq)
at = 100-gc									#AT/GC content of SARS sequence
print('GC content of SARS genome: ' + str(gc))
print('AT content of SARS genome: ' + str(at))

sequence = ''
count =  {'x<=20':0,'20<x<=40':0,'40<x<=60':0,'60<x<=80':0,'80<x<=100':0}  
beyond = 0

for i in range(len(covid[0].seq)):
	num = random.randint(1,100) 						#random sequence
	if num>=1 and num<30:
		sequence = sequence + 'G' 
	if num>=30 and num<60:
		sequence = sequence + 'C'
	if num>=60 and num<80:
		sequence = sequence + 'A'
	if num>=80 and num<=100:
		sequence = sequence + 'T'

sequence = Seq(sequence)	
print('GC content of random sequence: ' + str(GC(sequence)))		 	#GC content of random sequence
print('AT content of random sequence: ' + str(100-GC(sequence)))

proteins = sequence.translate()  						#translate
proteins = proteins.split('*')
 
for i in range(len(proteins)):
	if len(proteins[i]) <= 20:
		count['x<=20'] = count['x<=20']+1
	elif len(proteins[i]) <= 40:
		count['20<x<=40'] = count['20<x<=40']+1
	elif len(proteins[i]) <= 60:						#cutoff
		count['40<x<=60'] = count['40<x<=60']+1
	elif len(proteins[i]) <= 80:
		count['60<x<=80'] = count['60<x<=80']+1
	elif len(proteins[i]) <= 100:
		count['80<x<=100'] = count['80<x<=100']+1
	else:
		beyond = beyond +1
print('x = length of amino acid sequence')
print(count)
print('x>100: '+str(beyond))	
