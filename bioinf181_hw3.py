import random
import copy
import itertools

class GraphAlgorithms:

    def openFile(self, filename):
        with open(filename) as f:
            content = f.readlines()
        content = [x.strip() for x in content]

        return content

    def splitArrows(self, openedFile):
        paths = {}
        for lines in openedFile:
            partline = lines.split(" -> ")
            splitends = partline[1].split(",")
            for ends in splitends:
                if partline[0] not in paths:
                    paths[partline[0]] = [ends]
                else:
                    paths[partline[0]].append(ends)
        return paths

    def generateKmerComposition(self, string, kmer):
        outputKmers = []
        length = len(string) - kmer + 1
        for i in range(length):
            outputKmers.append(string[i:i+kmer])
        return outputKmers

    def reconstructString(self, sequences):
        seqLength = len(sequences[0])
        masterSequence = sequences.pop(0)
        while len(sequences) > 0:
            for seq in sequences:
                if masterSequence[len(masterSequence)-seqLength+1: len(masterSequence)] == seq[0:seqLength-1]:
                    masterSequence += seq[seqLength-1:seqLength]
                    sequences.remove(seq)
        return masterSequence

    def overlap(self, sequences):
        seqLen = len(sequences[0])
        overlapMap = {}
        for seq in sequences:
            for seq2 in sequences:
                if seq != seq2:
                    if seq[1:seqLen] == seq2[0:seqLen-1]:
                        overlapMap[seq] = seq2
        return overlapMap

    def generateDeBrujin(self, sequence, kmer):
        sequences = self.generateKmerComposition(sequence,kmer-1)
        overlapMap = {}
        for seq in sequences:
            for seq2 in sequences:
                if seq[1:kmer-1] == seq2[0:kmer-2]:
                    if seq not in overlapMap:
                        overlapMap[seq] = [seq2]
                    elif seq2 not in overlapMap[seq]:
                        overlapMap[seq].append(seq2)
        # for a in overlapMap:
        #     print a, "->",
        #     for b in overlapMap[a]:
        #         if len(overlapMap[a]) == 1:
        #             print b
        #         else:
        #             if b != overlapMap[a][-1]:
        #                 print b+",",
        #             else:
        #                 print b
        return overlapMap

    def generateDeBrujinWithSequences(self, sequences):
        kmerComps = {}
        for seqs in sequences:
            if seqs[0:len(seqs)-1] not in kmerComps:
                kmerComps[seqs[0:len(seqs)-1]] = [seqs[1:len(seqs)]]
            else:
                kmerComps[seqs[0:len(seqs) - 1]].append(seqs[1:len(seqs)])

        # for a in kmerComps:
        #     print a, "->",
        #     for b in range(len(kmerComps[a])):
        #         if len(kmerComps[a]) == 1:
        #             print kmerComps[a][b]
        #         else:
        #             if b != len(kmerComps[a]) - 1:
        #                 print kmerComps[a][b] + ",",
        #             else:
        #                 print kmerComps[a][b]
        return kmerComps

    def eulerianCycle(self, graph):
        masterSequences = []
        curSeq = []
        copyGraph = copy.deepcopy(graph)
        iniNode = random.choice(copyGraph.keys())
        curSeq.append(iniNode)
        cont = True

        while cont == True:
            nextNode = random.choice(copyGraph[iniNode])
            curSeq.append(nextNode)
            copyGraph[iniNode].remove(nextNode)
            if len(copyGraph[iniNode]) == 0:
                del copyGraph[iniNode]
            iniNode = nextNode
            if iniNode not in copyGraph:
                masterSequences.append(curSeq)
                if len(copyGraph) == 0:
                    cont = False
                else:
                    iniNode = random.choice(copyGraph.keys())
                    curSeq = []
                    curSeq.append(iniNode)

        while len(masterSequences) > 1:
            matchingNum = ''
            curList = ''
            for k in masterSequences[0]:
                for j in range(1,len(masterSequences),1):
                    if k in masterSequences[j]:
                        matchingNum = k
                        curList = j
            # print masterSequences, curList, matchingNum
            rotatedString = self.rotateSequence(masterSequences[curList],matchingNum)
            self.insertSequence(masterSequences[0], rotatedString)
            masterSequences.pop(curList)

        # retstring = ""
        # for u in range(len(masterSequences[0])):
        #     if u != len(masterSequences[0]) - 1:
        #         retstring += masterSequences[0][u]
        #         retstring += '->'
        #     else:
        #         retstring += masterSequences[0][u]
        # print retstring
        return masterSequences[0]

    def rotateSequence(self, seq, num):
        newseq = seq
        while newseq[0] != num:
            newseq.pop(0)
            newseq.append(seq[0])
        return newseq

    def insertSequence(self, seq1, seq2):
        newstring = seq1
        startindex = seq1.index(seq2[0]) + 1
        seq2.pop(0)
        for a in seq2:
            seq1.insert(startindex, a)
            startindex +=1
        return newstring

    def eulerianPath(self, graph):
        keys = set()
        inouts = {}
        for j in graph:
            keys.add(j)
            for k in graph[j]:
                keys.add(k)
        for key in keys:
            inouts[key] = 0
        for m in graph:
            for n in graph[m]:
                inouts[n] += 1
                inouts[m] -= 1
        for vals in inouts:
            if inouts[vals] == -1:
                startnode = vals
            if inouts[vals] == 1:
                endnode = vals

        if endnode not in graph:
            graph[endnode] = [startnode]
        else:
            graph[endnode].append(startnode)
        newseq = self.eulerianCycle(graph)
        rotatedseq = self.modRotateSequence(newseq, startnode, endnode)
        rotatedseq.pop(len(rotatedseq)-1)

        # retstring = ""
        # for u in range(len(rotatedseq)):
        #     if u != len(rotatedseq) - 1:
        #         retstring += rotatedseq[u]
        #         retstring += '->'
        #     else:
        #         retstring += rotatedseq[u]
        # print retstring
        return rotatedseq

    def modRotateSequence(self, seq, startnode, endnode):
        newseq = seq
        while (newseq[len(seq)-1] != startnode) or (newseq[len(seq)-2] != endnode):
            newseq.pop(0)
            newseq.append(seq[0])
        return newseq

    def stringReconstruction(self, patterns):
        debruj = self.generateDeBrujinWithSequences(patterns)
        path = self.eulerianPath(debruj)
        return self.reconstructString(path)

    def circularString(self, k):
        stringpermutes = range(0,2**k,1)
        for strings in range(len(stringpermutes)):
            binstring = bin(stringpermutes[strings])
            binstring = binstring[2:len(binstring)]
            while len(binstring) < k:
                binstring = '0' + binstring
            stringpermutes[strings] = binstring

        debrujin = self.generateDeBrujinWithSequences(stringpermutes)
        eulerpath = self.eulerianCycle(debrujin)
        pathstring = ""
        for nums in eulerpath:
            pathstring += nums[len(nums)-1]
        pathstring = pathstring[0:len(pathstring)-1]
        return pathstring


    def splitBar(self, openedFile):
        pairs = []
        for lines in openedFile:
            partline = lines.split("|")
            pairs.append(partline)
        return pairs

    def stringReconstructionReadPairs(self, readpairs, k, d):
        kmers = k
        while len(readpairs) > 1:
            i = 0
            while i < len(readpairs):
                suf1 = readpairs[i][0]
                suf2 = readpairs[i][1]
                j = 0
                while j < len(readpairs):
                    pref1 = readpairs[j][0]
                    pref2 = readpairs[j][1]
                    if suf1[len(suf1)-kmers+1:len(suf1)] == pref1[0:kmers-1] and suf2[len(suf2)-kmers+1:len(suf2)] == pref2[0:kmers-1]:
                        readpairs[i][0] += pref1[kmers-1:]
                        readpairs[i][1] += pref2[kmers-1:]
                        readpairs.remove([pref1,pref2])
                    j+=1
                i+=1
        finalsequence = readpairs[0][0] + readpairs[0][1][-(kmers+d):]
        return finalsequence

    def generateContigs(self, patterns):
        contigdb= self.generateDeBrujinWithSequences(patterns)
        keys = set()
        inouts = {}
        for j in contigdb:
            keys.add(j)
            for k in contigdb[j]:
                keys.add(k)
        for key in keys:
            inouts[key] = 0
        for m in contigdb:
            for n in contigdb[m]:
                inouts[n] += 1
                inouts[m] -= 1
        i = 0
        while i < len(patterns):
            pat1 = patterns[1]
            u = 0
            while u < len(patterns):
                pat2 = patterns[u]
                print pat1, pat2
                if pat2 != patterns[i]:
                    if pat2[0:len(pat2) - 1] == pat1[1:]:
                        # if inouts[pat1[0:len(pat2) - 1]] == 0:
                            patterns[i] += pat2[len(pat2) - 1:]
                            patterns.remove(pat2)
                            i = 0
                u += 1
            i += 1
        print patterns
        print inouts

    def gappedGenome(self, patterns, k, d):
        firstPats = []
        secondPats = []
        for pats in patterns:
            firstPats.append(pats[0])
            secondPats.append(pats[1])
        prefixString = self.newReconstructString(firstPats)
        # print prefixString
        suffixString = self.newReconstructString(secondPats)
        # print suffixString
        # for i in range(k+d+1, len(prefixString), 1):
            # if prefixString[i] != suffixString[i-k-d]:
            #     return None
        finalString = prefixString + suffixString[len(suffixString)-(k+d):]
        return finalString

    def newReconstructString(self, patterns):
        masterstring = patterns[0]
        # print masterstring
        patterns.pop(0)
        while len(patterns) > 0:
            masterstring += patterns[0][len(patterns[0])-1:]
            patterns.pop(0)
        return masterstring


g = GraphAlgorithms()
f = g.openFile("input")

# print g.eulerianCycle(g.splitArrows(f))
# g.eulerianPath(g.splitArrows(f))

# sequences = ["CTTA","ACCA","TACC","GGCT","GCTT","TTAC"]
# print g.stringReconstruction(f)

# print g.circularString(8)

# g.generateContigs(f)

# print g.gappedGenome(g.splitBar(f), 50, 200)

print g.stringReconstructionReadPairs(g.splitBar(f), 30, 100)


# retstrings = g.generateKmerComposition("", 50)
# retstrings.sort()
# for s in retstrings:
#     print s

# print g.reconstructString(f)

# sequences = ["ATGCG", "GCATG", "CATGC", "AGGCA", "GGCAT"]
# overlaps = g.overlap(f)
# for a in overlaps:
#     print a, "->", overlaps[a]

# sequence = "AAGATTCTCTAC"
# print g.generateDeBrujin(sequence, 4)

# sequences = ["GAGG","CAGG","GGGG","GGGA","CAGG","AGGG","GGAG"]
# g.generateDeBrujinWithSequences(f)

# sequences = ["CTTA","ACCA","TACC","GGCT","GCTT","TTAC"]
# g.stringReconstruction(sequences)


