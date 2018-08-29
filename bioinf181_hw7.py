class CombinatorialPatternMatching:

    def openFile(self, filename):
        with open(filename) as f:
            content = f.readlines()
        content = [x.strip() for x in content]
        return content

    def trieConstruction(self, seqs):
        trie = {}
        nodes = [0]
        paths = []

        for seq in seqs:
            curNode = 0
            if curNode not in trie:
                trie[curNode] = []
            for i in range(len(seq)):
                curSymbol = seq[i]
                inCurNode = False
                nextNode = 0
                if curNode not in trie:
                    trie[curNode] = []
                for a in trie[curNode]:
                    if nodes[a] == curSymbol:
                        inCurNode = True
                        nextNode = a
                if inCurNode == True:
                    curNode = nextNode
                else:
                    nodes.append(seq[i])
                    trie[curNode].append(len(nodes)-1)
                    paths.append((curNode, len(nodes)-1, seq[i]))
                    curNode = len(nodes)-1
        # Output
        # print paths
        # for b in paths:
        #     printString = ""
        #     printString += str(b[0]) + "->" + str(b[1]) + ":" + str(b[2])
        #     print printString
        return trie, nodes

    def modTrieConstruction(self, seqs, numArray):
        trie = {}
        nodes = [0]
        paths = []

        for seq in range(len(seqs)):
            curNode = 0
            if curNode not in trie:
                trie[curNode] = []
            for i in range(len(seqs[seq])):
                curSymbol = seqs[seq][i]
                inCurNode = False
                nextNode = 0
                if curNode not in trie:
                    trie[curNode] = []
                for a in trie[curNode]:
                    if nodes[a] == curSymbol:
                        inCurNode = True
                        nextNode = a
                if inCurNode == True:
                    curNode = nextNode
                else:
                    nodes.append(seqs[seq][i])
                    trie[curNode].append(len(nodes)-1)
                    paths.append((curNode, len(nodes)-1, seqs[seq][i]))
                    curNode = len(nodes)-1
            if curNode not in trie:
                trie[curNode] = []
            nodes.append(numArray[seq])
            trie[curNode].append(len(nodes) - 1)

        print trie, nodes
        return trie, nodes

    def trieMatching(self, text, patterns):
        textArray = []
        newtext = text
        curNum = 0
        numArray = []
        for a in range(len(newtext)):
            textArray.append(newtext)
            newtext = newtext[1:]
            numArray.append(curNum)
            curNum += 1

        trie, nodes = self.modTrieConstruction(textArray, numArray)
        # print trie, nodes, patterns
        # print nodes[161]
        output = []
        for seq in range(len(patterns)):
            curNode = 0
            # print patterns[seq]
            for i in range(len(patterns[seq]) + 1):
                curSymbol = patterns[seq][i-1]
                inCurNode = False
                nextNode = 0
                for a in trie[curNode]:
                    if nodes[a] == curSymbol:
                        inCurNode = True
                        nextNode = a

                if inCurNode == True:
                    curNode = nextNode

            # print curNode, trie[curNode]
            # for z in trie[curNode]:
            #     nxt = z
            #     while
            #     if
            # print curNode, nodes[curNode]

    def burrowsWheelerTransform(self, seq):
        rearrangements = []
        bwt = []
        curSeq = seq
        for i in range(len(seq)):
            lastSymb = curSeq[-1:]
            curSeq = lastSymb + curSeq[:-1]
            rearrangements.append(curSeq)
        sortedRearr = sorted(rearrangements)
        for a in sortedRearr:
            lastSym = a[-1:]
            bwt.append(lastSym)

        # output
        outline = ""
        for b in bwt:
            outline += b
        print outline

        return bwt

    def triMatch(self, sequences, patterns):
        entList = set()
        for pats in patterns:
            for i in range(len(sequences) - len(pats) + 1):
                if pats == sequences[i:i+len(pats)]:
                    entList.add(i)

        entList = sorted(entList)
        for a in entList:
            print a,

    def suffixArray(self, seq):
        rearrangements = []
        suffixArr = []
        curSeq = seq
        # suffixes = []
        for i in range(len(seq)):
            lastSymb = curSeq[-1:]
            curSeq = lastSymb + curSeq[:-1]
            rearrangements.append(curSeq)
        sortedRearr = sorted(rearrangements)

        for a in sortedRearr:
            indexPos = 0
            for u in range(len(a)):
                if a[u] == "$":
                    indexPos = u
            suffixArr.append(len(a) - indexPos - 1)
            # suffixes.append(a[0:indexPos+1])

        # output
        # outstring = ""
        # for b in suffixArr:
        #     outstring += str(b) + ", "
        # print outstring
        return suffixArr

    def multiplePatternMatchingSA(self, text, pattern, suffixArray):
        # for z in range(len(suffixArray)):
        #     print z, suffixArray[z], text[suffixArray[z]:]
        minIndex = 0
        maxIndex = len(text)
        while minIndex < maxIndex:
            midIndex = (minIndex + maxIndex)/2
            if pattern > text[suffixArray[midIndex]:suffixArray[midIndex] + len(pattern)]:
                minIndex = midIndex + 1
            else:
                maxIndex = midIndex
        first = minIndex
        maxIndex = len(text)

        while minIndex < maxIndex:
            midIndex = ((minIndex + maxIndex)/2)
            # print minIndex, maxIndex
            if pattern < text[suffixArray[midIndex]:suffixArray[midIndex] + len(pattern)]:
                maxIndex = midIndex
            else:
                minIndex = midIndex + 1
        last = maxIndex
        return first, last
        # print text[suffixArray[first + 1]:]

    # def suffixTreeCombine(self, trie, nodes):
    #     stack = [0]
    #     while len(stack) > 0:
    #         curNode = stack.pop()
    #         if curNode in trie:
    #             if len(trie[curNode]) == 1:
    #                 if (isinstance(nodes[trie[curNode][0]], int) == False):
    #                     nodes[curNode] = nodes[curNode] + nodes[trie[curNode][0]]
    #                     print nodes
    #                     for z in trie[nodes[trie[curNode][0]]]:
    #                         trie[curNode].append(trie[z])
    #                     trie[curNode].remove(0)
    #             for a in trie[curNode]:
    #                 stack.append(a)
g = CombinatorialPatternMatching()
f = g.openFile("input")

# pats = []
# for a in range(1, len(f)):
#     pats.append(f[a])
# g.trieMatching(f[0], pats)

# g.burrowsWheelerTransform(f[0])

# pats = []
# for a in range(1, len(f)):
#     pats.append(f[a])
# g.triMatch(f[0], pats)

# g.suffixArray(f[0])

print f
instring = f[0]
instring += "$"
suffixArray = g.suffixArray(instring)
positions = set()

for a in range(1, len(f)):
    first, last = g.multiplePatternMatchingSA(instring, f[a], suffixArray)
    for b in range(first,last):
        positions.add(suffixArray[b])
positions = sorted(positions)
for c in positions:
    print str(c),


# suffixArray, suffixes = g.suffixArray(f[0])
# trie, nodes = g.modTrieConstruction(suffixes, suffixArray)
# g.suffixTreeCombine(trie, nodes)