word = "ALABAMA"
word = word + "$"
anagrams  = []

for i in range(0, len(word), 1):
    char = word[0]
    word = word[1 :] + char
    anagrams.append(word)

anagrams.sort()

code = ""

for key in anagrams:
   code += key[ len(key) - 1]

print(code) 
