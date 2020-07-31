import string
from Bio.Seq import Seq
from Bio import SeqIO


filepath1 = "/share/SARS/betacov.fasta"
#Beta Coronavirus
filepath2 = "/share/SARS/SARS-urbani.fasta"
#SARS Coronavirus from 2003
filepath3 = "/share/SARS/mers.fasta"
#MERS Coronavirus
filepath4 = "/share/SARS/hku-bat.fasta"
# Bat1 coronavirus
filepath5 = "/share/SARS/nature-2013.fasta"
#Bat2 Coronavirus
filepath6 = "/share/SARS/pangolin-cov.fasta"
#Coronaviruspangolin
filepath7 = "/share/SARS/SARS-2020.fasta"       
#SARS-2Wuhan Coronavirus from Wuhan
filepath8 = "/share/SARS/covid-patient-CA.fasta"
#SARS2CAL1 corona virus from CA first patient
filepath9 = "/share/SARS/hku5-bat.fasta"
#Bat3
filepath10 ="/share/SARS/hedgehog-cov.fasta"
#hedgehog
def read_frames(filepath: str):
	reads = list(SeqIO.parse(filepath, "fasta"))
	for i in reads:
		rf1 = i.seq[0:]
		rf2 = i.seq[1:]
		rf3 = i.seq[2:]

	return rf1, rf2, rf3

def  protein_reader(s, Cut_off_count: int):
	p = s.transcribe()
	p = p.translate()
	counter = 0
	res = 0
	for i in range(0,len(p)) :
		if p[i] != '*':
			counter = counter + 1
		elif counter >= Cut_off_count and p[i] == '*':
			res = res + 1
			startLocation = i - counter
			endLocation = i-1
			print(startLocation, endLocation, counter, p[startLocation: endLocation])
			counter = 0
		else:
			counter = 0
	if counter >= Cut_off_count :
		res = res + 1
		startLocation = len - 1 - counter
		endLocation = len
		print(startLocation, endLocation, counter,  p[startLocation: endLocation])
	return res

cut_off = 100
r1, r2, r3 = read_frames(filepath1)
c1 = protein_reader(r1, cut_off)
c2 = protein_reader(r2, cut_off)
c3 = protein_reader(r3, cut_off)
print("Beta Corona Virus", c1, c2, c3, c1+c2+c3)

r1, r2, r3 = read_frames(filepath2)
c1 = protein_reader(r1, cut_off)
c2 = protein_reader(r2, cut_off)
c3 = protein_reader(r3, cut_off)
print("SARS Coronavirus from 2003",c1, c2, c3, c1+c2+c3)

r1, r2, r3 = read_frames(filepath3)
c1 = protein_reader(r1, cut_off)
c2 = protein_reader(r2, cut_off)
c3 = protein_reader(r3, cut_off)
print("MERS Coronavirus",c1, c2, c3, c1+c2+c3)

r1, r2, r3 = read_frames(filepath4)
c1 = protein_reader(r1, cut_off)
c2 = protein_reader(r2, cut_off)
c3 = protein_reader(r3, cut_off)
print("Bat Coronavirus",c1, c2, c3, c1+c2+c3)

r1, r2, r3 = read_frames(filepath5)
c1 = protein_reader(r1, cut_off)
c2 = protein_reader(r2, cut_off)
c3 = protein_reader(r3, cut_off)
print("Bat Coronavirus 2",c1, c2, c3, c1+c2+c3)

r1, r2, r3 = read_frames(filepath6)
c1 = protein_reader(r1, cut_off)
c2 = protein_reader(r2, cut_off)
c3 = protein_reader(r3, cut_off)
print("Coronavirus from Pangolin", c1, c2, c3, c1+c2+c3)

r1, r2, r3 = read_frames(filepath7)
c1 = protein_reader(r1, cut_off)
c2 = protein_reader(r2, cut_off)
c3 = protein_reader(r3, cut_off)
print("SARS-2 Coronavirus from Wuhan", c1, c2, c3, c1+c2+c3)

r1, r2, r3 = read_frames(filepath8)
c1 = protein_reader(r1, cut_off)
c2 = protein_reader(r2, cut_off)
c3 = protein_reader(r3, cut_off)
print("SARS-2 Coronavirus from first patient CA", c1, c2, c3, c1+c2+c3)

r1, r2, r3 = read_frames(filepath9)
c1 = protein_reader(r1, cut_off)
c2 = protein_reader(r2, cut_off)
c3 = protein_reader(r3, cut_off)
print("Bat Corona Virus 3", c1, c2, c3, c1+c2+c3)

r1, r2, r3 = read_frames(filepath10)
c1 = protein_reader(r1, cut_off)
c2 = protein_reader(r2, cut_off)
c3 = protein_reader(r3, cut_off)
print("hedgehog", c1, c2, c3, c1+c2+c3)
