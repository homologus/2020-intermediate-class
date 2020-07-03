#!/usr/bin/env python3
# Smith Waterman Algorythm (DNA Allignment AI)
# Written by CoolCat467 07/02/2020

NAME = 'Smith Waterman Algorythm'
__version__ = '0.0.1'

class SMAlgorythm(object):
    def __init__(self, sequence1, sequence2):
        self.seqA = str(sequence1).upper()
        self.seqB = str(sequence2).upper()
        self.n = len(self.seqA)
        self.m = len(self.seqB)
        self.ntides = ['A', 'T', 'G', 'C']
        self.subMatrix = self.getAllignmentDictionary()
        self.S = lambda a, b: self.subMatrix[(a, b)]
        # Linear Width
        self.w1 = 2
        self.W = lambda k: k*self.w1
        # Affine Width
        self.v = 5#Opening gap penalty
        self.u = 1#Gap Extention penalty
##        self.W = lambda k: (self.u*k) + self.v if self.v > 0 and self.u > 0 else 999999
        self.scoreMatrix = [[0 for i in range(self.m+1)] for i in range(self.n+1)]
    
    def __repr__(self):
        return "SMAlgorythm()"
    
    def getAllignmentDictionary(self):
        shift = lambda x: x[-1:] + x[:-1]
        x = list(self.ntides)
        tides = list(x)
        for i in range(len(self.ntides)-1):
            x = shift(x)
            tides += x
        isSame = lambda x, y: 3 if x==y else -3
        return {t:isSame(*t) for t in zip(tides, self.ntides*len(self.ntides))}
    
    def score(self):
        a = self.seqA
        b = self.seqB
        H = list(self.scoreMatrix)
        for k in range(self.n):
            for l in range(self.m):
                H[k][0] = 0
                H[0][l] = 0
        # main
        for i in range(1, self.n):
            for j in range(1, self.m):
                H[i][j] = max(H[i-1][j-1] + self.S(a[i], b[j]),#Alligning Score
                              H[i-k][j] - self.W(k),#Score if A is at the end of a gap of length k
                              H[i][j-l] - self.W(l),#Score if B is at the end of a gap of length l
                              0)
        self.scoreMatrix = H
        self.k = k
        self.l = l
    
    def getStartIdx(self):
        cmax = 0
        idx = None
        for i in range(self.n+1):
            for ii in range(self.m+1):
                if self.scoreMatrix[i][ii] > cmax:
                    cmax = self.scoreMatrix[i][ii]
                    idx = [i, ii]
        self.startIdx = tuple(idx)
    
    def getPath(self):
        i, j = tuple(self.startIdx)
        path = [tuple(self.startIdx)]
        while self.scoreMatrix[i][j] != 0:
            frm = (self.scoreMatrix[i-1][j-1] + self.S(self.seqA[i], self.seqB[j]),#Alligning Score
                   self.scoreMatrix[i-self.k][j] - self.W(self.k),#Score if A is at the end of a gap of length k
                   self.scoreMatrix[i][j-self.l] - self.W(self.l),#Score if B is at the end of a gap of length l
                   0)
            idxs = {frm[0]:(i-1, j-1),
                    frm[1]:(i-self.k, j),
                    frm[2]:(i, j-self.l)}
            m = max(frm)
            if m != 0:
                i, j = idxs[m]
                path.append(idxs[m])
                continue
            break
        self.path = [[self.seqA[i], self.seqB[j]] for i, j in reversed(path)]
    
    @classmethod
    def pathToString(cls, path):
        relative = []
        for i, ii in path:
            relative.append([i, '|' if i == ii else ' ', ii])
        data = ['', '', '']
        for r in range(len(relative)):
            i, rel, ii = relative[r]
            data[0] += i
            data[1] += rel
            data[2] += ii
        return ''.join([''.join(data[i])+'\n' for i in range(3)])[:-1]
    
    def allign(self, toString=False):
        self.score()
        self.getStartIdx()
        self.getPath()
        if toString:
            return self.__class__.pathToString(self.path)
        return self.path
    pass

if __name__ == '__main__':
    seq1 = 'catatgcgcattatgat'
##    seq1 = seq1 + ''.join([i for i in reversed(seq1)])
    seq2 = 'gatatgatcattatgct'
    print(seq1)
    print(seq2)
    print(SMAlgorythm(seq1, seq2).allign(True))
