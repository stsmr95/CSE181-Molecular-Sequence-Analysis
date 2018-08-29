from string import maketrans

class CombinatorialAlgorithms:

    def openFile(self, filename):
        with open(filename) as f:
            content = f.readlines()
        content = [x.strip() for x in content]
        return content

    def greedySorting(self, array):
        index = 0
        while index < len(array):
            if array[index] != index + 1:
                u = 0
                while abs(array[u]) != index + 1:
                    u+=1
                array = self.reverseString(array,index,u)
                self.printArray(array)
            else:
                index += 1

    def numBreakpoints(self, array):
        numbp = 0
        if array[0] != 1:
            numbp += 1
        curNum = array[0]
        index = 0
        while index < len(array) - 1:
            if curNum >= 0:
                if array[index + 1] != curNum + 1:
                    numbp += 1
                curNum = array[index+1]
                index += 1
            else:
                if array[index + 1] != curNum + 1:
                    numbp += 1
                curNum = array[index+1]
                index += 1
        print numbp

    def printArray(self,array):
        retString = "("
        for a in array:
            if a > 0:
                retString += "+" + str(a) + " "
            else:
                retString += str(a) + " "
        retString = retString[:-1]
        retString += ")"
        print retString,

    def printArray2(self,array):
        retString = "("
        for a in array:
            retString += str(a) + ", "
        retString = retString[:-2]
        retString += "),"
        print retString,

    def reverseString(self, array, pos1, pos2):
        newArray = []
        index = 0
        while index < pos1:
            newArray.append(array[index])
            index += 1
        backwardsindex = pos2
        while backwardsindex >= pos1:
            newArray.append(-1 * array[backwardsindex])
            backwardsindex -= 1
            index += 1
        while index < len(array):
            newArray.append(array[index])
            index += 1
        return newArray

    def parseStringToArray(self, string):
        parsedString = string[0].split("(")[1].split(")")[0].split()
        for a in range(len(parsedString)):
            parsedString[a] = int(parsedString[a])
        return parsedString

    def parseStringToArrays(self, string):
        s = string[0].split("(")
        s.pop(0)
        for a in range(len(s)):
            s[a] = s[a].split(")")[0].split()
            for b in range(len(s[a])):
                s[a][b] = int(s[a][b])
        return s

    def parseStringToGraph(self, string):
        s = string[0].split("(")
        s.pop(0)
        for a in range(len(s)):
            s[a] = s[a].split(")")[0].split(", ")
            for b in range(len(s[a])):
                s[a][b] = int(s[a][b])
        return s

    def sharedKmers(self, kmers, string1, string2):
        dic = {}
        for a in range(len(string1) - kmers + 1):
            if string1[a:a+kmers] not in dic:
                dic[string1[a:a+kmers]] = [a]
            else:
                newEntry = dic[string1[a:a+kmers]]
                newEntry.append(a)
                dic[string1[a:a + kmers]] = newEntry

            # if

        for b in range(len(string2) - kmers + 1):
            # print string2[b:b+kmers]
            revComp = self.reverseComplement(string2[b:b+kmers])
            if string2[b:b+kmers] in dic:
                for c in dic[string2[b:b+kmers]]:
                    outstr = "("
                    outstr += str(c) + ", " + str(b) + ")"
                    print outstr

            if revComp in dic:
                for d in dic[revComp]:
                    outstr = "("
                    outstr += str(d) + ", " + str(b) + ")"
                    print outstr

    def reverseComplement(self, sequence):
        intab = "ACGT"
        outtab = "TGCA"
        trantab = maketrans(intab, outtab)

        str = sequence.translate(trantab)
        return str[::-1]

    def chromosomeToCycle(self,chromo):
        cycle = []
        for i in range(0, len(chromo)):
            if chromo[i] > 0:
                cycle.append(2 * (chromo[i]) - 1)
                cycle.append(2 * (chromo[i]))
            else:
                cycle.append(-2 * (chromo[i]))
                cycle.append(-2 * (chromo[i]) - 1)
        return cycle

    def cycleToChromosome(self,cycle):
        chromo = []
        for j in range(0,len(cycle)/2):
            if cycle[(2*j)] < cycle[(2*j)+1]:
                chromo.append((cycle[(2*j)]+1)/2)
            else:
                chromo.append(-1*(cycle[(2*j)])/2)
        return chromo

    def coloredEdges(self, chromos):
        edges = []
        for a in range(len(chromos)):
            cycle = self.chromosomeToCycle(chromos[a])
            cycle.append(cycle[0])
            for i in range(1, len(chromos[a])+1):
                tuple = (cycle[2*i-1], cycle[2*i])
                edges.append(tuple)
        return edges

    def tupleOrigin(self, array):
        ori = []
        if array[0]%2 == 1:
            ori.append(-1*(array[0]+1)/2)
        else:
            ori.append((array[0] + 1) / 2)

        if array[1]%2 == 0:
            ori.append(-1*((array[1])/2))
        else:
            ori.append((array[1]+1)/2)

        return ori

    def graphToGenome(self, graph):
        oriGraph = []
        for a in range(len(graph)):
            oriGraph.append(self.tupleOrigin(graph[a]))
        allChromos = []
        while len(oriGraph) > 0:
            currentChromo = []
            tempArray = oriGraph.pop(0)
            for z in tempArray:
                currentChromo.append(z)
            while currentChromo[len(currentChromo)-1] != currentChromo[0]:
                for b in oriGraph:
                    if b[0] == currentChromo[len(currentChromo)-1]:
                        for d in b[1:]:
                            currentChromo.append(d)
                        oriGraph.remove(b)
            retChromo = currentChromo[:-1]
            allChromos.append(retChromo)

        # for y in allChromos:
        #     self.printArray(y)

        return allChromos
        # print oriGraph
        # for z in range(len(oriGraph)):
        #     oriGraph[z] = oriGraph[z][1:]
        # print oriGraph
        # print graph

    def twoBreakGenomeGraph(self, graph, i, j, a, b):
        firstappend = []
        secondappend = []
        if [i,j] in graph:
            graph.remove([i,j])
            firstappend = ([i,a])
        elif [j,i] in graph:
            graph.remove([j,i])
            firstappend = [a,i]

        if [a,b] in graph:
            graph.remove([a,b])
            secondappend = [j,b]
        elif [b,a] in graph:
            graph.remove([b,a])
            secondappend = [b,j]

        graph.append(firstappend)
        graph.append(secondappend)
        return graph

    def twoBreakGenomeTuple(self, parsedString, splitter):
        retArr = g.tupleToArray(g.coloredEdges(parsedString))
        arr = self.twoBreakGenomeGraph(retArr, int(splitter[0]), int(splitter[1]), int(splitter[2]), int(splitter[3]))
        retGenome = self.graphToGenome(self.arrayToTuple(arr))
        for a in retGenome:
            self.printArray(a)
        return retGenome

    def tupleToArray(self, tuples):
        retArray = []
        for a in tuples:
            retArray.append([a[0], a[1]])
        return retArray

    def arrayToTuple(self, arr):
        retTupArray = []
        for a in arr:
            retTupArray.append((a[0], a[1]))
        return retTupArray

g = CombinatorialAlgorithms()
f = g.openFile("input")
"""Problems 50-51
parsedString = g.parseStringToArray(f)

# g.greedySorting(parsedString)
g.numBreakpoints(parsedString)
"""

# Problem 54
# g.sharedKmers(int(f[0]), f[1], f[2])

# Problem 55
# parsedString = g.parseStringToArray(f)
# e = g.chromosomeToCycle(parsedString)
# for a in e:
#     print a,

# Problem 56
# parsedString = g.parseStringToArray(f)

# g.printArray(g.cycleToChromosome(parsedString))

# parsedString = g.parseStringToArrays(f)
# print g.coloredEdges(parsedString)

# parsedGraph = g.parseStringToGraph(f)
# print g.graphToGenome(parsedGraph)

# parsedGraph = g.parseStringToGraph(f)
# splitter = f[1].split(", ")
# arr = g.twoBreakGenomeGraph(parsedGraph, int(splitter[0]),int(splitter[1]),int(splitter[2]), int(splitter[3]))
# for a in arr:
#     g.printArray2(a)

# print f
parsedString = g.parseStringToArrays(f)
splitter = f[1].split(", ")

g.twoBreakGenomeTuple(parsedString, splitter)
# print parsedString


