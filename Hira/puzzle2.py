
word = "DVAAN$E"
#word = "EHLWA$"
#word = "ALM$BAAA"


l1 = []
for i in word:
    l1.append(i)

l2 = l1.copy()
l2.sort()

answer = l2[l1.index("$")]
l2.pop(l1.index("$"))
l1.pop(l1.index("$"))

for i in range(0, len(word) - 1, 1):
    txt = answer[i]
    l1Index = l1.index(txt)
    char =  l2[l1Index]
    if char == "$" and len(answer) < len(word) - 2 :
#         char = l2[l1.index(txt, l2.index("$"))] 
         indices = [i for i, x in enumerate(l1) if x == "A"]
         l1Index = indices[1]
    char =  l2[l1Index]
    answer = answer + char
    l2.pop(l1Index)
    l1.pop(l1Index)

answer = answer[ 0 : len(answer) - 1 ] 
print(answer)


