f = open("/share/Ecoli/GCA_000005845.2_ASM584v2_genomic.fna")
lines = f.readlines()
genome = ""
for i in range(1, len(lines)):
	for j in range(len(lines[i]) - 1):
		genome = genome + lines[i][j]
for i in range(0 , 3):
	extra = len(genome) - i % 11
	fragments = ""
	fragcount = 0
	fraglength = 0
	for j in range(i, len(genome) - extra, 11):
		site = genome[j:j+10]
		
