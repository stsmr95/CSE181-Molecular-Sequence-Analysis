class DynamicProgramming2:
    global blosmatrix
    blosmatrix = [[0 for a in range(25)] for b in range(25)]

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

    def topologicalOrder(self, list):
        # find incoming edges
        orderedList = []
        incEdges = {}
        for m in list.keys():
            if m not in incEdges:
                incEdges[m] = 0
            for n in list[m]:
                if n not in incEdges:
                    incEdges[n] = 1
                else:
                    incEdges[n] += 1

        complete = False
        while not complete:
            complete = True
            for c in incEdges:
                if incEdges[c] != 0:
                    complete = False

            for a in incEdges:
                if a not in orderedList:
                    if incEdges[a] == 0:
                        orderedList.append(a)
                        if a in list:
                            for b in list[a]:
                                incEdges[b] -= 1
        printstr = ""
        for z in orderedList:
            printstr = printstr + z + ", "
        printstr = printstr[:-2]
        print printstr
        return orderedList


    def multipleSequenceAlignment(self, seq1, seq2, seq3):
        n = len(seq1)
        m = len(seq2)
        l = len(seq3)
        pLMatrix = [[[0 for a in range(l + 1)] for b in range(m + 1)] for z in range(n+1)]
        bTMatrix = [[[0 for c in range(l + 1)] for d in range(m + 1)] for e in range(n+1)]
        pLMatrix[0][0][0] = 0
        # for u in range(1, m + 1):
        #     pLMatrix[0][u][0] = pLMatrix[0][u - 1][0] + 1
        # for v in range(1, n + 1):
        #     pLMatrix[v][0][0] = pLMatrix[v - 1][0][0] +1
        # for w in range(1, l + 1):
        #     pLMatrix[0][0][w] = pLMatrix[0][0][w - 1] +1
        # pLMatrix[0][0][0] = 0
        # print pLMatrix
        for x in range(1, n + 1):
            bTMatrix[x][0][0] = 1
        for y in range(1, m + 1):
            bTMatrix[0][y][0] = 2
        for z in range(1, l + 1):
            bTMatrix[0][0][z] = 3
        for a in range(1, n+1):
            for b in range(1, m+1):
                bTMatrix[a][b][0] = 4
        for c in range(1, n+1):
            for d in range(1, l+1):
                bTMatrix[c][0][d] = 5
        for e in range(1, m+1):
            for f in range(1, l+1):
                bTMatrix[0][e][f] = 6
        bTMatrix[0][0][0] = 0
        for i in range(1, n + 1, 1):
            for j in range(1, m + 1, 1):
                for k in range(1, l+1, 1):
                    # print seq1[i-1], seq2[j-1], seq3[k-1]
                    diag = pLMatrix[i - 1][j - 1][k - 1]
                    if seq1[i-1] == seq2[j-1] == seq3[k-1]:
                        diag = pLMatrix[i - 1][j - 1][k - 1] + 1
                    pLMatrix[i][j][k] = max(pLMatrix[i - 1][j][k],
                                        pLMatrix[i][j - 1][k],
                                        pLMatrix[i][j][k-1],
                                        pLMatrix[i - 1][j - 1][k],
                                        pLMatrix[i - 1][j][k - 1],
                                        pLMatrix[i][j - 1][k - 1],
                                        diag)

                    if pLMatrix[i][j][k] == pLMatrix[i - 1][j][k]:
                        bTMatrix[i][j][k] = 1
                    elif pLMatrix[i][j][k] == pLMatrix[i][j - 1][k]:
                        bTMatrix[i][j][k] = 2
                    elif pLMatrix[i][j][k] == pLMatrix[i][j][k - 1]:
                        bTMatrix[i][j][k] = 3
                    elif pLMatrix[i][j][k] == pLMatrix[i - 1][j - 1][k]:
                        bTMatrix[i][j][k] = 4
                    elif pLMatrix[i][j][k] == pLMatrix[i - 1][j][k - 1]:
                        bTMatrix[i][j][k] = 5
                    elif pLMatrix[i][j][k] == pLMatrix[i][j - 1][k - 1]:
                        bTMatrix[i][j][k] = 6
                    elif pLMatrix[i][j][k] == diag:
                        bTMatrix[i][j][k] = 7

        # print pLMatrix
        # for a in bTMatrix:
        #     print a
        # print pLMatrix
        return bTMatrix, seq1, seq2, seq3, n, m, l, pLMatrix[n][m][l]


    def outputMSAP(self, bTMatrix, seq1, seq2, seq3, n, m, l):
        # print bTMatrix
        iniString = ""
        secString = ""
        terString = ""
        totalLength = 0

        while bTMatrix[n][m][l] != 0 :
            if bTMatrix[n][m][l] == 1:
                n -= 1
                iniString = seq1[n] + iniString
                secString = '-' + secString
                terString = '-' + terString
            elif bTMatrix[n][m][l] == 2:
                m -= 1
                iniString = '-' + iniString
                secString = seq2[m] + secString
                terString = '-' + terString
            elif bTMatrix[n][m][l] == 3:
                l -= 1
                iniString = '-' + iniString
                secString = '-' + secString
                terString = seq3[l] + terString
            elif bTMatrix[n][m][l] == 4:
                n -= 1
                m -= 1
                iniString = seq1[n] + iniString
                secString = seq2[m] + secString
                terString = '-' + terString
            elif bTMatrix[n][m][l] == 5:
                n -= 1
                l -= 1
                iniString = seq1[n] + iniString
                secString = '-' + secString
                terString = seq3[l] + terString
            elif bTMatrix[n][m][l] == 6:
                m -= 1
                l -= 1
                iniString = '-' + iniString
                secString = seq2[m] + secString
                terString = seq3[l] + terString
            else:
                n -= 1
                m -= 1
                l -= 1
                iniString = seq2[m] + iniString
                secString = seq1[n] + secString
                terString = seq3[l] + terString
            totalLength += 1

        # while n == 0 and m != 0 and l != 0:
        #     iniString

        # print n, m ,l, totalLength
        print iniString
        print secString
        print terString
        return iniString, secString, terString

    def blosum62(self, sym1, sym2):
        return blosmatrix[ord(sym1)-65][ord(sym2)-65]


    def setUpBlosum62(self):
        openedFile = g.openFile("BLOSUM62")
        columns = openedFile[0].split()
        for a in range(1,len(openedFile)):
            line = openedFile[a].split()
            row = ord(line[0])-65
            for b in range(0, len(columns)):
                column = ord(columns[b])-65
                for c in range(1,len(line)):
                    blosmatrix[row][column] = int(line[b+1])


    def affineGapAlignment(self, seq1, seq2):
        n = len(seq1)
        m = len(seq2)

        midMatrix = [[0 for a in range(m + 1)] for b in range(n + 1)]
        bTMatrixl = [[0 for c in range(m + 1)] for d in range(n + 1)]
        bTMatrixu = [[0 for c in range(m + 1)] for d in range(n + 1)]
        bTMatrixm = [[0 for c in range(m + 1)] for d in range(n + 1)]

        upperMatrix = [[0 for a in range(m + 1)] for b in range(n + 1)]
        lowerMatrix = [[0 for a in range(m + 1)] for b in range(n + 1)]

        for i in range(1, n + 1, 1):
            for j in range(1, m + 1, 1):

                lowerMax = [lowerMatrix[i - 1][j] - 1, midMatrix[i - 1][j] - 11]
                lowerMatrix[i][j] = max(lowerMax)
                bTMatrixl[i][j] = lowerMax.index(lowerMatrix[i][j])

                upperMax = [upperMatrix[i][j-1] - 1, midMatrix[i][j-1] - 11]
                upperMatrix[i][j] = max(upperMax)
                bTMatrixu[i][j] = upperMax.index(upperMatrix[i][j])

                midMax = [midMatrix[i - 1][j - 1] + self.blosum62(seq1[i - 1], seq2[j - 1]),
                            lowerMatrix[i][j],
                            upperMatrix[i][j]]
                midMatrix[i][j] = max(midMax)
                bTMatrixm[i][j] = midMax.index(midMatrix[i][j])

        # for a in bTMatrixm:
        #    print a
        # for a in bTMatrix:
        #     print a
        # print pLMatrix
        return bTMatrixl, bTMatrixm, bTMatrixu, seq1, seq2, n, m, midMatrix[n][m]

    def outputAGAP(self, bTMatrixl, bTMatrixm, bTMatrixu, seq1, seq2, n, m):
        iniString = ""
        secString = ""
        currentMatrix = 1
        # print n, m, bTMatrixm[n][m]

        while n != 0 and m != 0:
            # lower level
            # print n, m, bTMatrixm[n][m], currentMatrix, iniString, secString
            if currentMatrix == 0:
                if bTMatrixl[n][m] == 0:
                    iniString = seq1[n-1] + iniString
                    secString = '-' + secString
                    n -= 1
                else:
                    n -= 1
                    iniString = seq1[n] + iniString
                    secString = '-' + secString
                    currentMatrix = 1

            # middle level
            elif currentMatrix == 1:
                if bTMatrixm[n][m] == 0:
                    iniString = seq1[n-1] + iniString
                    secString = seq2[m-1] + secString
                    n -= 1
                    m -= 1
                elif bTMatrixm[n][m] == 1:
                    currentMatrix = 0
                else:
                    currentMatrix = 2

            # upper level
            else:
                if bTMatrixu[n][m] == 0:
                    iniString = '-' + iniString
                    secString = seq2[m-1] + secString
                    m -= 1
                else:
                    m -= 1
                    iniString = '-' + iniString
                    secString = seq2[m] + secString
                    currentMatrix = 1

        print iniString
        print secString

        return iniString, secString

    def midEdges(self, seq1, seq2):
        n = len(seq1)
        m = len(seq2)
        mid = m/2
        mu = 5
        column1 = [0 for a in range(n + 1)]
        column2 = [0 for b in range(n + 1)]
        bTcolumn1 = [0 for c in range(n + 1)]
        bTcolumn2 = [0 for d in range(n + 1)]

        for j in range(1, mid + 1, 1):
            for i in range(1, n + 1, 1):
                diag = column1[i - 1] + self.blosum62(seq1[i - 1], seq2[j - 1])
                column2[i] = max(column2[i-1] - mu,
                                     column1[i] - mu,
                                     diag)

                if column2[i] == column1[i] - mu:
                    bTcolumn1[i] = 3
                elif column2[i] == column2[i-1] - mu:
                    bTcolumn1[i] = 1
                elif column2[i] == diag:
                    bTcolumn1[i] = 2

            column1 = column2
            column2 = [0 for j in range(n + 1)]

        column3 = [0 for e in range(n + 1)]
        column4 = [0 for f in range(n + 1)]


        for j in range(m-1, mid - 1, -1):
            for i in range(n-1, -1, -1):
                diag2 = column1[i + 1] + self.blosum62(seq1[i], seq2[j])
                column3[i] = max(column3[i + 1] - mu,
                                 column4[i] - mu,
                                 diag2)
                if column3[i] == column4[i] - mu:
                    bTcolumn2[i] = 1
                elif column3[i] == column3[i + 1] - mu:
                    bTcolumn2[i] = 3
                elif column3[i] == diag2:
                    bTcolumn2[i] = 2
            column4 = column3
            column3 = [0 for h in range(n + 1)]

        newcolumn = [0 for e in range(n + 1)]
        for f in range(n+1):
            newcolumn[f] = column1[f] + column4[f]
        maxVal = max(newcolumn)
        # print "MIDNODE:", newcolumn.index(maxVal), mid
        print newcolumn

        curRow = newcolumn.index(maxVal)
        while bTcolumn2[curRow] == 3:
            curRow += 1
        if bTcolumn2[curRow] == 2:
            curRow += 1

        retString = "("
        retString = retString + str(newcolumn.index(maxVal)) + ", " + str(mid) + ") "
        retString = retString + "(" + str(curRow) + ", " + str(mid+1) + ")"
        print retString

        return newcolumn.index(maxVal), mid


g = DynamicProgramming2()
f = g.openFile("input")
# g.topologicalOrder(g.splitArrows(f))

# bTMatrix, seq1, seq2, seq3, n, m, l, maxi = g.multipleSequenceAlignment(f[0], f[1], f[2])
# print maxi
# g.outputMSAP(bTMatrix, seq1,seq2,seq3,n,m,l)

g.setUpBlosum62()
# bTMatrixl, bTMatrixm, bTMatrixu, seq1, seq2, n, m, maxi = g.affineGapAlignment(f[0],f[1])
# print maxi
# g.outputAGAP(bTMatrixl, bTMatrixm, bTMatrixu,seq1,seq2,n,m)

g.midEdges(f[0], f[1])