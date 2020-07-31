f = open("words.fasta")
lines = f.readlines()
seqList = []
for i in lines:
	for j in range(0, len(i), 17):
		seqList.append(i[j:j+15])
seq = seqList[0]
temp = []
for i in range(1, len(seqList)):
	for j in range(4, 15):
		if(seqList[i][0:j] == seq[len(seq)-j:len(seq)]):
			seq = seq + seqList[i][j:15]
	for j in range(1, 12):
		if(seqList[i][j:15] == seq[0:15-j]):
			seq = seqList[i][0:j] + seq
		else:
			temp.append(seqList[i])
for i in range(len(temp)):
        for j in range(4, 15):
                if(temp[i][0:j] == seq[len(seq)-j:len(seq)]):
                        seq = seq + temp[i][j:15]
        for j in range(1, 12):
                if(temp[i][j:15] == seq[0:15-j]):
                        seq = temp[i][0:j] + seq
print(seq)
