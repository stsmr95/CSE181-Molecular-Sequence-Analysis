from bioinf181_hw1 import PatternCompute
import random
import copy

class MotifSearch:

    def weightedChoice(self, weights):
        totals = []
        runningTotal = 0

        for w in weights:
            runningTotal += w
            totals.append(runningTotal)

        rnd = random.random() * runningTotal
        for i, total in enumerate(totals):
            if rnd < total:
                return i

    def motifEnumeration(self, sequences, kmer, mismatch):
        patterns = set()
        pc = PatternCompute()

        allpats = set()
        u = 0
        while u < 4**kmer:
            allpats.add(pc.numberToPattern(u,kmer))
            u += 1

        for pat in allpats:
            allcontains = True

            for seq in sequences:
                returnedlist = pc.approximatePattern(pat,seq,mismatch)
                if len(returnedlist) == 0:
                    allcontains = False

            if allcontains == True:
                patterns.add(pat)

        return patterns


    def medianString(self, sequences, kmer):
        distance = 9999
        median = ""
        pc = PatternCompute()

        allpats = set()
        u = 0
        while u < 4**kmer:
            allpats.add(pc.numberToPattern(u,kmer))
            u += 1

        for pats in allpats:
            totalHamDistance = 0
            for seq in sequences:
                i = 0
                minHamDistance = 9999
                while i < len(seq) - kmer:
                    hamdis = pc.hammingDistance(pats, seq[i: i+kmer])
                    if hamdis < minHamDistance:
                        minHamDistance = hamdis
                    i += 1
                totalHamDistance += minHamDistance
            if totalHamDistance < distance:
                distance = totalHamDistance
                median = pats
        return median

    def distanceBetweenPatternAndStrings(self, pattern, sequences):
        pc = PatternCompute()
        totalHamDistance = 0
        for seq in sequences:
            i = 0
            minHamDistance = 9999
            while i < len(seq) - len(pattern):
                hamdis = pc.hammingDistance(pattern, seq[i: i + len(pattern)])
                if hamdis < minHamDistance:
                    minHamDistance = hamdis
                i += 1
            totalHamDistance += minHamDistance
        return totalHamDistance


    def profileMostProbableKmer(self, profile, sequence, kmer):
        highestprob = float(0.0)
        hpseq = sequence[0:kmer]
        i = 0
        while i <= len(sequence) - kmer:
            u = 0
            totprob = float(1.0)
            while u < kmer:
                prob = float(0.0)
                if sequence[u + i] == "A":
                    prob = profile[0][u]
                elif sequence[u + i] == "C":
                    prob = profile[1][u]
                elif sequence[u + i] == "G":
                    prob = profile[2][u]
                else:
                    prob = profile[3][u]
                totprob = float(totprob * prob)
                u += 1
            if totprob > highestprob:
                highestprob = totprob
                hpseq = sequence[i:i+kmer]
            i += 1
        return hpseq

    def updateCountChart(self, pattern, chart):
        x = 0
        while x < len(pattern):
            if pattern[x] == "A":
                chart[0][x] += 1
            elif pattern[x] == "C":
                chart[1][x] += 1
            elif pattern[x] == "G":
                chart[2][x] += 1
            else:
                chart[3][x] += 1
            x += 1
        return chart

    def createProfile(self, chart):
        profile = [[0.0 for a in range(len(chart[0]))] for b in range(4)]
        count = chart[0][0] + chart[1][0] + chart[2][0] + chart[3][0]
        y = 0
        while y < len(chart[0]):
            profile[0][y] = float(chart[0][y]) / count
            profile[1][y] = float(chart[1][y]) / count
            profile[2][y] = float(chart[2][y]) / count
            profile[3][y] = float(chart[3][y]) / count
            y += 1

        return profile

    def createProfileWithPseudocounts(self, chart):
        profile = [[0.0 for a in range(len(chart[0]))] for b in range(4)]
        for nucleotide in range(len(chart)):
            for letter in range(len(chart[nucleotide])):
                chart[nucleotide][letter] += 1

        count = chart[0][0] + chart[1][0] + chart[2][0] + chart[3][0]
        y = 0
        while y < len(chart[0]):
            profile[0][y] = float(chart[0][y]) / count
            profile[1][y] = float(chart[1][y]) / count
            profile[2][y] = float(chart[2][y]) / count
            profile[3][y] = float(chart[3][y]) / count
            y += 1

        return profile

    def score(self, motifs):
        scoreArray = []
        numMotifs = len(motifs)
        i = 0
        while i < len(motifs[0]):
            arrayScore = 0
            numA = 0
            numC = 0
            numG = 0
            numT = 0
            u = 0
            while u < numMotifs:
                if motifs[u][i] == "A":
                    numA += 1
                elif motifs[u][i] == "C":
                    numC += 1
                elif motifs[u][i] == "G":
                    numG += 1
                else:
                    numT += 1
                u+=1
            arrayScore = numA + numC + numG + numT - max(numA, numC, numG, numT)
            # if consensus[i] == "A":
            #     arrayScore = numC + numG + numT
            # elif consensus[i] == "C"
            #     arrayScore = numA + numG + numT
            # elif consensus[i] == "G"
            #     arrayScore = numA + numC + numT
            # else
            #     arrayScore = numA + numC + numG
            scoreArray.append(arrayScore)
            i += 1
        totalscore = 0
        for sc in scoreArray:
            totalscore += sc
        return totalscore

    def randomizedMotifSearch(self, sequences, kmer):
        ms = MotifSearch()
        motifs = []
        bestmotifs = []
        for u in range(len(sequences)):
            rannum = random.randrange(0, len(sequences[0])-kmer+1, 1)
            appseq = sequences[u][rannum:rannum+kmer]
            motifs.append(appseq)
            bestmotifs.append(appseq)

        countMatrix = [[0 for a in range(kmer)] for b in range(4)]
        while 1 == 1:
            for line in motifs:
                countMatrix = ms.updateCountChart(line, countMatrix)
            profile = ms.createProfileWithPseudocounts(countMatrix)
            for o in range(len(sequences)):
                motifs[o] = ms.profileMostProbableKmer(profile, sequences[o], kmer)
            if ms.score(motifs ) < ms.score(bestmotifs):
                bestmotifs = motifs
            else:
                return bestmotifs

    def greedyMotifSearch(self, sequences, kmer):
        ms = MotifSearch()
        bestmotifs = []
        for strings in sequences:
            bestmotifs.append(strings[0:kmer])
        i = 0
        while i < len(sequences[0]) - kmer:
            countMatrix = [[0 for a in range(kmer)] for b in range(4)]
            motif0 = sequences[0][i:i+kmer]
            motifs = [None for c in range(len(sequences))]
            motifs[0] = motif0
            countMatrix = ms.updateCountChart(motif0,countMatrix)
            for u in range(1,len(sequences),1):
                profile = ms.createProfile(countMatrix)
                newmotif = ms.profileMostProbableKmer(profile,sequences[u],kmer)
                motifs[u] = newmotif
                countMatrix = ms.updateCountChart(newmotif, countMatrix)
            if ms.score(motifs) < ms.score(bestmotifs):
                bestmotifs = motifs
            i+=1
        return bestmotifs

    def gibbsSampler(self, sequences, kmer, N):
        ms = MotifSearch()
        bestmotifs = [None for a in range(len(sequences))]
        motifs = [None for a in range(len(sequences))]
        for sequence in range(len(sequences)):
            ranmotif = random.randrange(0, len(sequences[0]) - kmer + 1, 1)
            bestmotifs[sequence] = sequences[sequence][ranmotif:ranmotif+kmer]
            motifs[sequence] = sequences[sequence][ranmotif:ranmotif + kmer]
        countMatrix = [[0 for a in range(kmer)] for b in range(4)]
        for j in range(1, N, 1):

            selectedkmernum = random.randrange(0, len(sequences), 1)
            for kmers in range(len(motifs)):
                if kmers != selectedkmernum:
                    countMatrix = ms.updateCountChart(motifs[kmers], countMatrix)
            profile = ms.createProfileWithPseudocounts(countMatrix)
            probabilities = [0 for a in range(len(sequences[selectedkmernum]) - kmer)]
            possibleoptions = [None for a in range(len(sequences[selectedkmernum]) - kmer)]
            for k in range(0, len(sequences[selectedkmernum]) - kmer):
                possibleoptions[k] = sequences[selectedkmernum][k:k+kmer]
                prob = 1.0
                for l in range(0,kmer,1):
                    if sequences[selectedkmernum][l+k] == "A":
                        prob *= profile[0][l]
                    elif sequences[selectedkmernum][l+k] == "C":
                        prob *= profile[1][l]
                    elif sequences[selectedkmernum][l+k] == "G":
                        prob *= profile[2][l]
                    else:
                        prob *= profile[3][l]
                probabilities[k] = prob
            totalprob = sum(probabilities)
            for m in range(len(probabilities)):
                probabilities[m] = float(probabilities[m] / (totalprob))
            choice = ms.weightedChoice(probabilities)
            motifs[selectedkmernum] = possibleoptions[choice]
            if ms.score(motifs) < ms.score(bestmotifs):
                bestmotifs = copy.deepcopy(motifs)
        return bestmotifs



m = MotifSearch()
"""sequences = {"CTTGTACTACCAGAGACGTTCAAGG", "GCGCGTGCCGCCATCCCCGGGGCAT", "CATCTATTGGATGACTCGCGCCCGG", "ATATACGTGGATTGACTCGACACAT", "TTTCGTCTATCACGGGACGAAATTG", "TTGGTCAGGGTGTACTCGTGACATT", "GACCACTGGGACTTGTTGTGGCGAC", "CGGGGGAGGTCCTATCACAAGGACT", "CAGGGAGGTGATTTGGTGGCGGAAA", "CTAGGCTGATTACATCATCTTTCAG"}
retpatterns = m.motifEnumeration(sequences, 5, 2)
for pats in retpatterns:
    print pats"""

"""sequences = {"CTGTGAGGGACGACTCACATGAACTGAATAATAGTGAGTACC","GAGCTCAAGCCGAATGAGAGGCAAGGGATGTAGTGGGAGCAA","AAGACGGTCCCGGGGAGGGGAGGCTCAAAGCAGTGGGCCTTT","GCCAGGGGGATGTAAAGTCCGCTTTTCAAGGAGAGATTCCCA","AGGGCGCGGACGGCAAAACCGACCTCGGGAGGGAAGCACAAC", "GTATGGATCCAGAAATGGTCCTTACTAGACGGGAAGTAATAT", "AATAGTTCCTTTATCCGACGTCGAGGGATGCTAGAACACACC", "GGGAAGGCGCATCCACAGAGGATGCTATGTCCGTCAAGGTCG", "ACCCGTCCGAGCCACTCACGGCCAAAGCCTAAAGGCGGGATG", "ATGTCTGGGATGCGCCTATGAATTCACCGCGGGAGACTGCCG"}
median = m.medianString(sequences, 6)
print median"""

"""profile = [[0.303, 0.364, 0.121, 0.212, 0.333, 0.182] , [0.182, 0.182, 0.242, 0.182, 0.091, 0.273], [0.333, 0.242, 0.303, 0.394, 0.364, 0.273], [0.182, 0.212, 0.333, 0.212, 0.212, 0.273]]
print m.profileMostProbableKmer(profile, "CAGACGTAGACTGCCGTCTCAGTCTGTTATGTCTCGTTGAACGCTGGGAATATGCTGATAAGGGGAGCAAAATCCGCAAGCGTGGTGGTCTCAGAGCATCGGTAGTTTTTACCTAGACGTATATGGCCCTCGACCGCAAGTATCGGCGGGACGTGTAGATCCAAAGTATCTACAAGTCTCTAGATTAGCATTGTTAGCGG", 6)
"""

"""profile = [[0.2, 0.2, 0.3, 0.2, 0.3],[0.4, 0.3, 0.1, 0.5, 0.1],[0.3, 0.3, 0.5, 0.2, 0.4],[0.1, 0.2, 0.1 ,0.1, 0.2]]
print m.profileMostProbableKmer(profile, "ACCTGTTTATTGCCTAAGTTCCGAACAAACCCAATATAGCCCGAGGGCCT", 5)"""

"""sequences = ["CTATCCCCACACCCGCCGATGACACGAAGCGCTTGTTGTGTTTTCATTGCAAGATCGCTAGGGAGACGTGACGTACCGTACGAGTGGGAGTGTATAGCCCCACGGCTGATGCGAGGTGATTTTAGGTTAGAATCAGATGTGCGCTCTCAGTAAGGT",
"GATGGGCGACGGGAGACTTCCGGTTTCAGTCCTTCCTGTTTTATCAGTGAGGGCGCTAATTGCACTATAACCTTCGTCAGCAGTTTCTAGGCATACGCAATGTACATTACCGCAGGGTGTTAAACAACCATTATCTTTATCCTTTCGAAGCCAGGT",
"TCGGGGGCCTTGTCGGCAAAGTTGAACTGGAACTAAATTCTGCCTTGTAATCAGGTGTAGAGTTTGTCCGCAAGCATCAAAGGTTTCTAATAACTCATCACGCGAACTTCAAAGGTAGGTAATGGAAGACACTATACGGGGTTGGACTAACTCTCC",
"GCTTATAACGTTATAAGCTTTCATTGTAGAAAAAGGAGATCATATTCCTGATAAGCGTCTCGCGTGTGGTGACGTCAGAATACCGGAGCGTCGAGATCTCAGCCAGGTTTCCTGGTGATCAGGTAAGTCCCGATATCTTCTAAAAATGAATTCACT",
"TCAATAGCGTTTCATCCCGTGGGATCTGAGCTAGGTCGAGGAAACGGCATGGTGTGTGCGTAGCCCCCTGAGTCCCGCAGTACCAGGCAGAAGACATAACAACGATAACACGGACCGAGGATATGAATGTCGAAGATCCCAGTAGACCTCAAGGTT"]


greed = m.greedyMotifSearch(sequences, 12)
for i in greed:
    print i"""

"""f = open('rosalind_ba2h.txt', 'r')
pattern = f.readline()

# print sequences
print m.distanceBetweenPatternAndStrings("GGCTT",sequences)"""

"""chart = [[1,1,2],[1,2,1],[1,0,0],[0,0,0]]"""

"""motifs = ["TCGGGGGTTTTT", "CCGGTGACTTAC", "ACGGGGATTTTC"]
print m.score(motifs)"""

# sequences = ["CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA", "GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG", "TAGTACCGAGACCGAAAGAAGTATACAGGCGT", "TAGATCAAGTTTCAGGTGCACGTCGGTGAACC", "AATCCACCAGCTCCACGTGCAATGTTGGCCTA"]

sequences = ["GTTGAATGTCGAGGACGGTGGCATCGTATCGTCTTAACAGTCGCCACTGATAACCCTTTTTTGGGTTTGAGCTGATCGGATAGGCCAAGGAACGAATAAACCGTTACATTATATTTGCCATTGCCCGGTATTCCGGGTGGTCGTGGTCTTCAGCAACCGTTGAATGTCGAGGA",
"CGGTGGCATCGTATCGTCTTAACAGTCGCCACTGATAACCCTTTTTTGGGTTTGAGCTGATCGGATAGGCCAAGGAACGAATAAACCGTTACATTATATGCTGAGAAATAAGCTTTGCCATTGCCCGGTATTCCGGGTGGTCGTGGTCTTCAGCAACCGTTGAATGTCGAGGA",
"GGATATCAATGGTGAGAGTAATGTTCGACGGTATCGTGAGCCAATGTCCGCAGATTCTCCGCACCATAATTTGTGCCAGCGCATAAGCTAAAAGCTACTCGCGTGATCCACATCGGCAAGCACCTTCATTTTCGGACCGGTCAGTCGGGCGAACAGAATAAAGCTACTGTAAA",
"AGTTTTTATCCACTGATAGCCAGGAAATAACGCACGGGCGCCACCAAGAAAAACAATCATGAACGGGCGTGAGAGTAGCGGAAAGCTGTGAAACGGGCAGGGCGTGCGTTAGTCCATGACTATAACAAACTCAGCCGAGTGTAATACCATGACCAGCCACCTTATTCCAGGGA",
"GGGCGGTATAAGCCCGCATATGTCAGACCACAGGGGTCTTTCGCGTTCTAACATATGGCTAACTGTCCCGAGTGTGTGAAGTAAGTCGACCATAACGGCACCAGACCTCAGGCCAGCAAGTTACTCAGACTGGGGTTAGAGGAAATAAGCACAAGAGTGCCTACCACATCTAG",
"AGGCTAAAGAATCTGAGGCGCGACGCGTAAGTCTGCGCGTATCTGCGAGCCGCTGTGTTTTGAAGTTTATGATATCCACGCCTCAAACATTCAGGAAATAAGCCGGGCTGTCTTGCTATATAGCTTGATTCCGGGTCGAGCTGCCGATAATGCGTACTTAAAGTAAATCGACC",
"AAACATATACGAGACACACTCATGTGTTAGAGCACAACGTTCAATTCCATTACAATATGATTATAGCTAACACTTAACTTCGGGCAAGGCTTGGCGAACGATAAGCAGAGGCTTTCAAGCTTTTTAGAAAGCACTGCCAGGAAAATGGCAGCTCACTTAGCTACACCAGGAGT",
"CATACTTTTCCTGAATGATTTATAGTTAAAATAGCTCAACTTGTCACGCCTTCCGTGTGGCCGCGTTATGTGTTATGGAGTGCCAATGAATAAGCTGTCCCGGACCATCCGACGGAGATTGGGTTACAAATTCGTGGGATGGGGAGTGGAGCCGAGTCGGGGAGTGAGACACA",
"ACCAGACTTGCTTGAACAGCTTCTATACTGTTGGGAAATAAGCCGCGTTCCGTTTCCGCCCTTCGTAAGATAGGAAACGTGGGCCGAGCTAATACAGCATCTCCTTAGAGTATAGATTCCTCCTACATAGGTCTCAGTTTACATGACCTGATCGAGAAGGTGATGCGGCCGCC",
"GGCTGCTTGATTCAGCACTCCAAAGAATATACGGCACCCCAGTCGACCCGCTCCTTAGCACCCATCCTCTCCAGCAACCTAACCTCATAATTTCTCTCACTATACTGGGTTAGTGCAGGGCCTAGCGGTGCCAATCAATAAGCCGTATTTTCCGAACATATAGGATTGCGGCT",
"TCGGAAGTAGTTGAGAATAGAAACGGCGAGCCATGGAACGTGGAATAGTCTCGACACTTCTACCTGTATCATCTTTTGAACAAAAACCTTTGCCCCTAAATAAGCAAGATATGCATTACATTAAGTAATGTAATTGAGAGCCTATACACCTCAGCAATCCCGACGCTCGGTCT",
"TAAAAGACTTTGTGTGAGGATGCCAGCCTATAAGCCAAAACAACAAGTCACCCCTACTCTACCTGTGCCCAACTGTACCCCCCTATGATGTACTCGCGGGATGCTGTCTCCCACGAGTCTAGGCGCGACAATTCCTCAAAGTCGACCGTCGTTAGTCGCACAGTCTAGCTGCA",
"CCTGCTCTCCGTCTTATGCCGGGATGCTGTGGTCCACCAACCTGATTAGCGCAGTGACTTCTCGTGCTCATAATAAGTGATATGACCATAAATAGTGTTACACGCTCGCTGAGAATGCAGTGAAATAAGCCGCATTGGGCGAATCCCCTGCGGGCTACACCTGCGCGGAGTCA",
"TGCCGATTGGCGGCGACGTGTGGTCCTTGTGCCAGGAAATGCCCACTCGCGTGGCAATCACAACTTGATAGTCGCCCCTTATAGCCAGATGAAACAATAACGGTTCTGGAAACGTTGCACATTAACGAACTATGTCCTCCGTCATAGGTACCAAACGCGGTGCCGAAGTCCGA",
"ATTGCACCAGCATATGACGGGACGTGACTTTAGCAATGGCGGTTCAGGGTGTATTTGGCATCTATGGGCTATACTTTAAGGGATATGCCAGGCTGTAAGCACTGCAGCGCGAGGAAAGCTCGGGATAATACCAGGCAGCACAAGTGGTGGTCGACGTATATATACCTATCCTG",
"AAGAAACACACAGGATGCGCATACCACCAGGAAATAAGACTTACTATCTGTCCGTACGCACGTCGTCGTAGCAGAGTGGTTCATGCGTGGCTGCAAATATTTTCCGTAAACTAGAAACGTACCATCGGCAGGACAAATGACTCCAGCAAGGGAAATCACTCATGATGGTTCCA",
"CGGGATTCGCAACGGTTATTCCATTGGCGGAAAAGCTTAGCCGGCAATGGGCAATCAGTTTTAGGCAACATCCCATGTAGGTCAGACCGTACCAGGAGGTGTGCTTAACTGGCAGATATTGCGGGTGGGCGCCCAATGCCAGGAAATATTGTGCGTTTTCGGTCTATGCACGA",
"GAGCATGCACGAATGAATCAGAAATACTTCATGGACTAACTTGTACCATATGGGGAATGCCGCCAAATAAGCGAAGATCATCTGGGCCATGTCGCTCCGGTAACGGAGAATCTGGCGGAAACGGTCGCAGTACGGGCGGACACATATTCCGGGACGGCCGCTAGGACTGTGGT",
"GCCTTAGGTATGCCAGGAAGCGAGCACACGCTACCTTCGAAAACAGTTTTCAAGTGAGAGCGACATCTTACGTCTCATGTGCTCCAGTGCCGGCTCCGTCCCGTATACAGTCCTACATAAGCTTACATCTCCGTCCTCACTCGTGTTGATAGTCACCGAGAATACCTGGGGAG",
"TTAGGAGAGCCCAACATAACTGAGCTAGCTCTCTGATTCCTTCCTCTCGGTATCACATCCTCAGACACATCTTATTACTTACCACTCCGATCTTGATAGAATTCTCCTTTCATTCCAGGCGGAGACCGAGTTAACGTCTGCCAGGATGCAAGCCTCATGCGTGTGTATAGGCG"]
bestmotifs = m.randomizedMotifSearch(sequences, 15)
bestscore = m.score(bestmotifs)
for u in range(2000):
    motiflist = m.randomizedMotifSearch(sequences,15)
    # print bestmotifs, bestscore
    if m.score(motiflist) < bestscore:
        bestmotifs = motiflist
        bestscore = m.score(motiflist)

print "Lowest score:", bestscore
for v in bestmotifs:
    print v


"""sequences = ["AATTGGCACATCATTATCGATAACGATTCGCCGCATTGCC", "GGTTAACATCGAATAACTGACACCTGCTCTGGCACCGCTC", "AATTGGCGGCGGTATAGCCAGATAGTGCCAATAATTTCCT", "GGTTAATGGTGAAGTGTGGGTTATGGGGAAAGGCAGACTG", "AATTGGACGGCAACTACGGTTACAACGCAGCAAGAATATT", "GGTTAACTGTTGTTGCTAACACCGTTAAGCGACGGCAACT", "AATTGGCCAACGTAGGCGCGGCTTGGCATCTCGGTGTGTG", "GGTTAAAAGGCGCATCTTACTCTTTTCGCTTTCAAAAAAA"]
sequences = ["CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA", "GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG", "TAGTACCGAGACCGAAAGAAGTATACAGGCGT", "TAGATCAAGTTTCAGGTGCACGTCGGTGAACC", "AATCCACCAGCTCCACGTGCAATGTTGGCCTA"]
numberMatrix = [[None] for b in range(len(sequences))]
for u in range(1000):
    motiflist = m.randomizedMotifSearch(sequences, 8)
    for i in range(len(motiflist)):
        numberMatrix[i].append(motiflist[i])

for v in range(len(numberMatrix)):
    print max(list(numberMatrix[v]), key=numberMatrix[v].count)"""

"""sequences = ["CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA",
"GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG",
"TAGTACCGAGACCGAAAGAAGTATACAGGCGT",
"TAGATCAAGTTTCAGGTGCACGTCGGTGAACC",
"AATCCACCAGCTCCACGTGCAATGTTGGCCTA"]"""

# sequences = []


"""numberMatrix = [[] for b in range(len(sequences))]
for u in range(20):
    motiflist = m.gibbsSampler(sequences,8, 200)
    for i in range(len(motiflist)):
        numberMatrix[i].append(motiflist[i])

for v in range(len(numberMatrix)):
    print max(list(numberMatrix[v]), key=numberMatrix[v].count)
"""

"""bestmotifs = m.gibbsSampler(sequences, 15, 2000)
bestscore = m.score(bestmotifs)
for u in range(20):
    motiflist = m.gibbsSampler(sequences,15, 2000)
    # print bestmotifs, bestscore
    if m.score(motiflist) < bestscore:
        bestmotifs = motiflist
        bestscore = m.score(motiflist)

for v in bestmotifs:
    print v
# print bestmotifs, bestscore"""