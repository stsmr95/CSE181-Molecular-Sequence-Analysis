class CombinatorialPatternMatching2:

    def openFile(self, filename):
        with open(filename) as f:
            content = f.readlines()
        content = [x.strip() for x in content]
        return content

    def inverseBWT(self, string):
        lastColumn = []
        for a in range(len(string)):
            lastColumn.append(string[a])

        firstColumn = sorted(lastColumn)
        lastCount = self.generateCountArray(lastColumn)
        firstCount = self.generateCountArray(firstColumn)
        oriString = " "
        curLastLetter = "$"
        for z in range(len(lastColumn)):
            if lastColumn[z] == curLastLetter:
                curLastPos = lastCount[z]
        nextLetter = ""
        nextPos = 0

        while oriString[len(oriString) - 1] != "$":
            for y in range(len(lastColumn)):
                if lastColumn[y] == curLastLetter:
                    if lastCount[y] == curLastPos:
                        nextLetter = firstColumn[y]
                        nextPos = firstCount[y]
            oriString += nextLetter
            curLastLetter = nextLetter
            curLastPos = nextPos

        return oriString[1:]

    def lastToFirst(self, column, place):
        lastCountArray = self.generateCountArray(column)
        firstArray = sorted(column)
        firstCountArray = self.generateCountArray(firstArray)
        letter = column[place]
        lastPlace = lastCountArray[place]
        for a in range(len(firstArray)):
            if firstArray[a] == letter:
                if firstCountArray[a] == lastPlace:
                    return a

    def generateCountArray(self, col):
        # firstColumn = sorted(lastColumn)
        countDict = {}
        countArray = []
        for a in range(len(col)):
            if col[a] not in countDict:
                countDict[col[a]] = 1
            else:
                countDict[col[a]] += 1
            countArray.append(countDict[col[a]])
        return countArray

    def multiplePatternMatching(self, text, pattern, suffixArray):
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

    def BWMatching(self, firstColumn, lastColumn, pattern):
        # print lastColumn
        # print firstColumn
        top = 0
        bottom = len(lastColumn) - 1
        # print top, bottom
        while top <= bottom:
            # print top, bottom
            if len(pattern) != 0:

                symbol = pattern[-1]
                pattern = pattern[:-1]
                contains = False
                lastPos = bottom
                firstPos = 9999999999
                # print symbol, top, bottom, pattern
                for a in range(top, bottom+1):
                    # print a, lastColumn[a]
                    if lastColumn[a] == symbol:
                        # print "!", a, lastColumn[a]
                        contains = True
                        lastPos = a
                        if a <= firstPos:
                            firstPos = a
                # print "FL",firstPos, lastPos
                if contains == True:
                    top = self.lastToFirst(lastColumn, firstPos)
                    bottom = self.lastToFirst(lastColumn, lastPos)
                    # print "New:", top, bottom
                else:
                    # print "cant find", pattern
                    return 0
            else:
                # print top, bottom
                return bottom - top + 1
        print top, bottom
        return 0

    def betterBWMatching(self, lastColumn, pattern):
        # print lastColumn
        # print firstColumn
        firstColumn = sorted(lastColumn)
        top = 0
        bottom = len(lastColumn) - 1
        # print top, bottom
        fOdict = self.firstOccurence(firstColumn)
        while top <= bottom:
            # print top, bottom
            if len(pattern) != 0:

                symbol = pattern[-1]
                pattern = pattern[:-1]
                contains = False
                for a in range(top, bottom + 1):
                    if lastColumn[a] == symbol:
                        contains = True
                if contains == True:
                    top = fOdict[symbol] + self.count(symbol, top, lastColumn)
                    bottom = fOdict[symbol] + self.count(symbol, bottom + 1, lastColumn) - 1
                else:
                    # print "cant find", pattern
                    return 0
            else:
                # print top, bottom
                return bottom - top + 1
        print top, bottom
        return 0

    def count(self,sym, num,column):
        count = 0
        for a in range(0, num):
            if column[a] == sym:
                count+=1
        return count

    def firstOccurence(self, column):
        dict = {}
        for a in range(len(column)):
            if column[a] not in dict:
                dict[column[a]] = a
        return dict

    def suffixArray(self, seq):
        rearrangements = []
        suffixArr = []
        curSeq = seq
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

        # output
        # outstring = ""
        # for b in suffixArr:
        #     outstring += str(b) + ", "
        # print outstring
        return suffixArr

g = CombinatorialPatternMatching2()
f = g.openFile("input")

# print g.inverseBWT(f[0])

# print g. lastToFirst(f[0], int(f[1]))
# lastColumn = []
# for a in range(len(f[0])):
#     lastColumn.append(f[0][a])
# print g.lastToFirst(lastColumn, int(f[1]))



lastColumn = []
for a in range(len(f[0])):
    lastColumn.append(f[0][a])
firstColumn = sorted(lastColumn)
allPats = f[1].split()
for b in allPats:
    print g.BWMatching(firstColumn, lastColumn, b),


# print g.firstOccurence("a", "$aaaaaabmnnnps")

#
# lastColumn = []
# for a in range(len(f[0])):
#     lastColumn.append(f[0][a])
# allPats = f[1].split()
# for b in allPats:
#     print g.betterBWMatching( lastColumn, b),
# #
# instring = f[0]
# instring += "$"
# suffixArray = g.suffixArray(instring)
# positions = set()
#
# for a in range(1, len(f)):
#     first, last = g.multiplePatternMatching(instring, f[a], suffixArray)
#     for b in range(first, last):
#         positions.add(suffixArray[b])
# positions = sorted(positions)
# for c in positions:
#     print str(c),