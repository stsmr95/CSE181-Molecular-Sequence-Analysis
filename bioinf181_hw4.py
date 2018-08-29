class DynamicProgramming:
    global blosmatrix
    blosmatrix = [[0 for a in range(25)] for b in range(25)]

    def openFile(self, filename):
        with open(filename) as f:
            content = f.readlines()
        content = [x.strip() for x in content]
        return content

    def parseCoins(self, openedFile):
        money = openedFile[0]
        coins = openedFile[1].split(',')
        for coin in range(len(coins)):
            coins[coin] = int(coins[coin])
        return int(money), coins


    def minimumNumCoins(self, money, coins):
        minCoinsArray = [0 for a in range(0, money+1)]
        currentMoney = 1
        while currentMoney < money+1:
            minCoinsArray[currentMoney] = 9999
            for change in coins:
                if currentMoney >= change:
                    if minCoinsArray[currentMoney-change] + 1 < minCoinsArray[currentMoney]:
                        minCoinsArray[currentMoney] = minCoinsArray[currentMoney-change] + 1
            currentMoney += 1
        print minCoinsArray[money]

    def parsePathArrays(self, openedFile):
        iniLine = openedFile[0]
        n = int(iniLine.split(" ")[0])
        m = int(iniLine.split(" ")[1])
        openedFile.pop(0)
        downMatrix = []
        rightMatrix = []
        for x in range(n):
            # print x, n
            array = openedFile[x].split(" ")
            for y in range(len(array)):
                array[y] = int(array[y])
            downMatrix.append(array)

        for a in range(n+1):
            # print a+1+n, m
            array = openedFile[a+1+n].split(" ")
            for b in range(len(array)):
                array[b] = int(array[b])
            rightMatrix.append(array)
        return n, m, downMatrix, rightMatrix

    def manhattanTourist(self, n, m, down, right):
        pLMatrix = [[0 for a in range(m+1)] for b in range(n+1)]
        for i in range(1,n+1,1):
            pLMatrix[i][0] = pLMatrix[i-1][0] + down[i-1][0]
        for j in range(1,m+1,1):
            pLMatrix[0][j] = pLMatrix[0][j-1] + right[0][j-1]
        for x in range(1, n+1, 1):
            for y in range(1,m+1,1):
                pLMatrix[x][y] = max(pLMatrix[x-1][y] + down[x-1][y], pLMatrix[x][y-1] + right[x][y-1])
        return pLMatrix

    def longestCommonSubsequence(self, openedFile):
        seq1 = openedFile[0]
        seq2 = openedFile[1]
        n = len(seq1)
        m = len(seq2)
        pLMatrix = [[0 for a in range(m + 1)] for b in range(n + 1)]
        bTMatrix = [[0 for c in range(m + 1)] for d in range(n + 1)]
        for i in range(1, n+1, 1):
            for j in range(1, m+1, 1):
                diag = pLMatrix[i-1][j-1]
                if seq1[i-1] == seq2[j-1]:
                    diag = pLMatrix[i-1][j-1] + 1
                pLMatrix[i][j] = max(pLMatrix[i-1][j], pLMatrix[i][j-1], diag)

                if pLMatrix[i][j] == pLMatrix[i-1][j]:
                    bTMatrix[i][j] = 3
                elif pLMatrix[i][j] == pLMatrix[i][j-1]:
                    bTMatrix[i][j] = 1
                elif pLMatrix[i][j] == diag:
                    bTMatrix[i][j] = 2
        print pLMatrix
        return bTMatrix, seq1, n, m

    def outputLCS(self, bTMatrix, seq1, n, m):
        curString = ""
        while n != 0 and m != 0:
            # print curString, n, m
            if bTMatrix[n][m] == 3:
                n -= 1
            elif bTMatrix[n][m] == 1:
                m -= 1
            else:
                n -= 1
                m -= 1
                curString = seq1[n] + curString
        print curString

    def parseDAG(self, openedFile):
        sourceNode = openedFile[0]
        sinkNode = openedFile[1]
        paths = []
        maxNodes = 0
        for a in range(2, len(openedFile)):
            curpath = openedFile[a].split("->")
            # print openedFile[a], curpath
            weightsplit = []
            weightsplit.append(curpath[1])
            weightsplit = weightsplit[0].split(":")
            curpath.pop(1)
            if maxNodes < int(weightsplit[0]):
                maxNodes = int(weightsplit[0])
            pathlist = curpath + weightsplit
            for b in range(len(pathlist)):
                pathlist[b] = int(pathlist[b])
            paths.append(pathlist)

        return int(sourceNode), int(sinkNode), paths, int(maxNodes)

    def longestDAGPath(self, sourceNode, sinkNode, paths, max):

        longestLengthArray = [0 for a in range(0, max+1)]
        prevNode = [0 for a in range(0, max+1)]
        totalWeights = [0 for a in range(0, max+1)]
        connectsToSource = [False for a in range(0, max+1)]
        connectsToSource[sourceNode] = True
        curNode = sourceNode
        while curNode != sinkNode:
            for a in range(len(paths)):
                if paths[a][0] == curNode:
                    # if paths[a][0] >= sourceNode:
                    #     if paths[a][1] < sinkNode:
                    if connectsToSource[paths[a][0]] == True:
                        # print paths[a], max
                        if longestLengthArray[paths[a][1]] < longestLengthArray[paths[a][0]] + 1:
                            longestLengthArray[paths[a][1]] = longestLengthArray[paths[a][0]] + 1
                            prevNode[paths[a][1]] = curNode
                            totalWeights[paths[a][1]] = totalWeights[paths[a][0]] + paths[a][2]
                            connectsToSource[paths[a][1]] = connectsToSource[paths[a][0]]
                    # paths.pop(a)
            curNode += 1

        backTrackNode = sinkNode
        bTArray = []
        while backTrackNode != sourceNode:
            bTArray.append(backTrackNode)
            backTrackNode = prevNode[backTrackNode]
        bTArray.append(sourceNode)
        bTArray.reverse()
        print totalWeights[sinkNode]

        retstring = ""
        for u in range(len(bTArray)):
            if u != len(bTArray) - 1:
                retstring += str(bTArray[u])
                retstring += '->'
            else:
                retstring += str(bTArray[u])
        print retstring

        return totalWeights[sinkNode], bTArray

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

    def setUpPAM250(self):
        openedFile = g.openFile("PAM250")
        columns = openedFile[0].split()
        for a in range(1,len(openedFile)):
            line = openedFile[a].split()
            row = ord(line[0])-65
            for b in range(0, len(columns)):
                column = ord(columns[b])-65
                for c in range(1,len(line)):
                    blosmatrix[row][column] = int(line[b+1])

    def globalAlignment(self, seq1, seq2):
        n = len(seq1)
        m = len(seq2)
        mu = 5
        pLMatrix = [[0 for a in range(m + 1)] for b in range(n + 1)]
        bTMatrix = [[0 for c in range(m + 1)] for d in range(n + 1)]
        pLMatrix[0][0] = 0
        for u in range(1, m + 1):
            pLMatrix[0][u] = pLMatrix[0][u - 1] - mu
        for v in range(1, n + 1):
            pLMatrix[v][0] = pLMatrix[v - 1][0] - mu
        pLMatrix[0][0] = 0
        for x in range(1, m + 1):
            bTMatrix[0][x] = 1
        for y in range(1, n + 1):
            bTMatrix[y][0] = 3
        bTMatrix[0][0] = 0
        for i in range(1, n + 1, 1):
            for j in range(1, m + 1, 1):
                diag = pLMatrix[i - 1][j - 1] + self.blosum62(seq1[i - 1], seq2[j - 1])
                pLMatrix[i][j] = max(pLMatrix[i - 1][j] - mu,
                                     pLMatrix[i][j - 1] - mu,
                                     diag)

                if pLMatrix[i][j] == pLMatrix[i - 1][j] - mu:
                    bTMatrix[i][j] = 3
                elif pLMatrix[i][j] == pLMatrix[i][j - 1] - mu:
                    bTMatrix[i][j] = 1
                elif pLMatrix[i][j] == diag:
                    bTMatrix[i][j] = 2

        # print pLMatrix
        # for a in bTMatrix:
        #     print a
        # print pLMatrix
        return bTMatrix, seq1, seq2, n, m, pLMatrix[n][m]


    def outputGAP(self, bTMatrix, seq1, seq2, n, m):
        # print bTMatrix
        curString = ""

        while bTMatrix[n][m] != 0:
            if bTMatrix[n][m] == 3:
                n -= 1
                curString = '-' + curString
            elif bTMatrix[n][m] == 1:
                m -= 1
                curString = seq2[m] + curString
            else:
                n -= 1
                m -= 1
                curString = seq2[m] + curString

        return curString

    def outputGAP2(self, bTMatrix, seq1, seq2, n, m):
        iniString = ""

        while bTMatrix[n][m] != 0:
            # print n,m, iniString
            if bTMatrix[n][m] == 3:
                n -= 1
                iniString = seq1[n] + iniString
            elif bTMatrix[n][m] == 1:
                m -= 1
                iniString = '-' + iniString
            else:
                n -= 1
                m -= 1
                iniString = seq1[n] + iniString

        return iniString


    def localAlignment(self, seq1, seq2):
        n = len(seq1)
        m = len(seq2)
        mu = 5
        pLMatrix = [[0 for a in range(m + 1)] for b in range(n + 1)]
        bTMatrix = [[0 for c in range(m + 1)] for d in range(n + 1)]
        pLMatrix[0][0] = 0
        maxi = 0
        o = p = 0

        # for u in range(1, m + 1):
        #     pLMatrix[0][u] = pLMatrix[0][u - 1] - mu
        # for v in range(1, n + 1):
        #     pLMatrix[v][0] = pLMatrix[v - 1][0] - mu
        # pLMatrix[0][0] = 0
        for x in range(1, m + 1):
            bTMatrix[0][x] = 1
        for y in range(1, n + 1):
            bTMatrix[y][0] = 3
        bTMatrix[0][0] = 0
        for i in range(1, n + 1, 1):
            for j in range(1, m + 1, 1):
                diag = pLMatrix[i - 1][j - 1] + self.blosum62(seq1[i - 1], seq2[j - 1])
                pLMatrix[i][j] = max(0,
                                     pLMatrix[i - 1][j] - mu,
                                     pLMatrix[i][j - 1] - mu,
                                     diag)

                if pLMatrix[i][j] > maxi:
                    o = i
                    p = j
                    maxi = pLMatrix[i][j]

                if pLMatrix[i][j] == pLMatrix[i - 1][j] - mu:
                    bTMatrix[i][j] = 3
                elif pLMatrix[i][j] == pLMatrix[i][j - 1] - mu:
                    bTMatrix[i][j] = 1
                elif pLMatrix[i][j] == diag:
                    bTMatrix[i][j] = 2

        # print pLMatrix
        # for a in pLMatrix:
        #     print a
        # print pLMatrix
        print maxi
        return bTMatrix, seq1, seq2, o, p, pLMatrix

    def outputLAP(self, bTMatrix, seq1, seq2, n, m, pLMatrix):
        # print bTMatrix
        curString = ""

        while pLMatrix[n][m] != 0:
            # if bTMatrix[n][m] != 0:
                if bTMatrix[n][m] == 3:
                    n -= 1
                    curString = '-' + curString
                elif bTMatrix[n][m] == 1:
                    m -= 1
                    curString = seq2[m] + curString
                else:
                    n -= 1
                    m -= 1
                    curString = seq2[m] + curString

        return curString

    def outputLAP2(self, bTMatrix, seq1, seq2, n, m, pLMatrix):
        iniString = ""

        while pLMatrix[n][m] != 0:
                # print n,m, iniString
                if bTMatrix[n][m] == 3:
                    n -= 1
                    iniString = seq1[n] + iniString
                elif bTMatrix[n][m] == 1:
                    m -= 1
                    iniString = '-' + iniString
                else:
                    n -= 1
                    m -= 1
                    iniString = seq1[n] + iniString

        return iniString

    def editDistance(self, string1, string2):
        disMatrix = [[0 for a in range(len(string2) + 1)] for b in range(len(string1) + 1)]
        for i in range(0, len(string2) + 1):
            disMatrix[0][i] = i
        for u in range(0, len(string1) + 1):
            disMatrix[u][0] = u
        for x in range(1, len(string1) +1):
            for y in range(1, len(string2) +1):
                if string1[x-1] == string2[y-1]:
                    disMatrix[x][y] = disMatrix[x-1][y-1]
                else:
                    # print x, y, y-1
                    insertDistance = disMatrix[x][y-1] + 1
                    deleteDistance = disMatrix[x-1][y] + 1
                    replaceDistance = disMatrix[x-1][y-1] + 1
                    disMatrix[x][y] = min(insertDistance, deleteDistance, replaceDistance)

        return disMatrix[len(string1)][len(string2)]

    def modGlobalAlignment(self, seq1, seq2):
        n = len(seq1)
        m = len(seq2)
        mu = 1
        pLMatrix = [[0 for a in range(m + 1)] for b in range(n + 1)]
        bTMatrix = [[0 for c in range(m + 1)] for d in range(n + 1)]
        pLMatrix[0][0] = 0
        for u in range(1, m + 1):
            pLMatrix[0][u] = pLMatrix[0][u - 1] - 1
        for v in range(1, n + 1):
            pLMatrix[v][0] = pLMatrix[v - 1][0] - 1
        pLMatrix[0][0] = 0
        for x in range(1, m + 1):
            bTMatrix[0][x] = 1
        for y in range(1, n + 1):
            bTMatrix[y][0] = 3
        bTMatrix[0][0] = 0
        for i in range(1, n + 1, 1):
            for j in range(1, m + 1, 1):
                diag = pLMatrix[i - 1][j - 1] - 1
                if seq1[i - 1] == seq2[j - 1]:
                    diag = pLMatrix[i - 1][j - 1] + 1
                pLMatrix[i][j] = max(pLMatrix[i - 1][j] - 1,
                                     pLMatrix[i][j - 1] - 1,
                                     diag)

                if pLMatrix[i][j] == pLMatrix[i - 1][j] - mu:
                    bTMatrix[i][j] = 3
                elif pLMatrix[i][j] == pLMatrix[i][j - 1] - mu:
                    bTMatrix[i][j] = 1
                elif pLMatrix[i][j] == diag:
                    bTMatrix[i][j] = 2

        # print pLMatrix
        # for a in bTMatrix:
        #     print a
        # print pLMatrix
        return bTMatrix, seq1, seq2, n, m, pLMatrix[n][m]

    def modLocalAlignment(self, seq1, seq2):
        n = len(seq1)
        m = len(seq2)
        mu = 1
        pLMatrix = [[0 for a in range(m + 1)] for b in range(n + 1)]
        bTMatrix = [[0 for c in range(m + 1)] for d in range(n + 1)]
        pLMatrix[0][0] = 0
        maxi = 0
        o = p = 0

        # for u in range(1, m + 1):
        #     pLMatrix[0][u] = pLMatrix[0][u - 1] - mu
        # for v in range(1, n + 1):
        #     pLMatrix[v][0] = pLMatrix[v - 1][0] - mu
        # pLMatrix[0][0] = 0
        for x in range(1, m + 1):
            bTMatrix[0][x] = 1
        for y in range(1, n + 1):
            bTMatrix[y][0] = 3
        bTMatrix[0][0] = 0
        for i in range(1, n + 1, 1):
            for j in range(1, m + 1, 1):
                diag = pLMatrix[i - 1][j - 1] - 1
                if seq1[i - 1] == seq2[j - 1]:
                    diag = pLMatrix[i - 1][j - 1] + 1
                pLMatrix[i][j] = max(0,
                                     pLMatrix[i - 1][j] - mu,
                                     pLMatrix[i][j - 1] - mu,
                                     diag)

                if pLMatrix[i][j] > maxi:
                    o = i
                    p = j
                    maxi = pLMatrix[i][j]

                if pLMatrix[i][j] == pLMatrix[i - 1][j] - mu:
                    bTMatrix[i][j] = 3
                elif pLMatrix[i][j] == pLMatrix[i][j - 1] - mu:
                    bTMatrix[i][j] = 1
                elif pLMatrix[i][j] == diag:
                    bTMatrix[i][j] = 2

        # print pLMatrix
        # for a in pLMatrix:
        #     print a
        # print pLMatrix
        # print maxi
        return bTMatrix, seq1, seq2, o, p, maxi, pLMatrix


    def fittingAlignmentProblem(self, string1, string2):
        print len(string1)
        print g.editDistance(string1, string2)
        highestscore = -999
        hin = him = 0
        readingframe = g.editDistance(string1, string2)
        histr1 = ""
        while readingframe > len(string2):
            for i in range(0, len(string1)-readingframe):
                newstring = string1[i:i+readingframe+1]
                # print "A"
                bTMatrix, seq1, seq2, n, m, maxi = self.modGlobalAlignment(newstring, string2)
                if maxi > highestscore:
                    highestscore = maxi
                    highestmatrix = bTMatrix
                    histr1 = newstring
                    hin = n
                    him = m
            readingframe -= 1

        print highestscore
        print self.outputGAP2(highestmatrix, histr1, string2, hin, him)
        print self.outputGAP(highestmatrix,histr1,string2,hin,him)




g = DynamicProgramming()
f = g.openFile("input")

# money, coins = g.parseCoins(f)
# g.minimumNumCoins(money, coins)

# n, m, down, right = g.parsePathArrays(f)
# manhat = g.manhattanTourist(n,m,down,right)
# print manhat[n][m]

# bTMatrix, seq1, n, m = g.longestCommonSubsequence(f)
# g.outputLCS(bTMatrix, seq1, n, m)

# source, sink, paths, max = g.parseDAG(f)
# g.longestDAGPath(source, sink, paths, max)


# print ord('A'), ord('C'), ord('Y')
# g.blosum62('A', 'B')

# g.setUpBlosum62()
# g.setUpPAM250()
# matrix, seq1, seq2, n, m, maxi = g.globalAlignment(f[0], f[1])
# matrix, seq1, seq2, n, m, pL= g.localAlignment(f[0], f[1])
# g.outputLAP2(matrix, seq1,seq2,n,m,pL)
# g.outputLAP(matrix, seq1,seq2,n,m,pL)

g.fittingAlignmentProblem(f[0], f[1])

# print maxi
# print maxi
# print seq1
# g.outputGAP2(matrix,seq1,seq2,n,m)
# g.outputGAP(matrix,seq1,seq2,n,m)


# print g.editDistance("PEGVKLHYNFSYFHEEPQVQLYCSGFDIRHVCWFFKDHREHLPGCCGKFAWFFRTTSREQNAKAASTMPPQDILGWFAMAMLRPEGSRQNMKMVYTHSQICQFDYDHFQEVIKQGKNDLSIPVCKKQELAIQCQHLTWKSVEKTGFENLWNTQPKPMESCDQYECRDMCCIRKKMIKRIFWHCWLCQHDFFATHWGAYILYTFPLDQLPKIEAHNHFQQKQPAQVNQHKLDRCFFRIAMVPGVREHIIFEPCRDCGCIAHFGSINVRTHNCCGRRGMPFADAGHTMQPKYTSEKVEWTYCYCREFNIPPPNWTMVATTYYDAICGVLCGMVCCFPHSVIYQNLSCWAEIEMRTPTVEEEHHVALSLGKEMIDICVMISTRQVTCMVFGICWNDCRRTHQDPSAIFTVDQQEGCDNPHCMEERWWEMTYMVENQYNTWKAKAVFACYSWHEGIELYVWDRMISSFPKFQPGDLWHAWNLKCNMDSESVSCPCYFPLPKANHCTRDMRVMATHTYWTIQQKWVPVQWKYGAHIWHEGCMNHVSPDKYFETESVGNIANITRDEHCWWVLVAHASFLPRFPPFPFREDIGKHKRRPITVHCAGWQVLCERGTIWRFLFHHDRIKTGNSSPLTAECQTVSMVCINCEAQYGFLSHFRSIGKISFQPNWWARSWLYIYLKYHTQVMNKCPDGAPNCHLSDDKYESMVMMWNLTYHASVHTYGMVYQKQITPFLSCLEVEHIVYRKGEFQCDYLLPTEQFYQCLSVGPMDGHVIPNPDLTLHFHHQTTFRRECCDPIYDILKGDFHYVDEHYFHLSAKFIQEICEESLQGFVSLPGYLAVNEVPFERAIDVDCRPHCKKMCLAWPCDIIHKRMPHFNTDKYEWKEQQLNPWMFQLNWKMPDRLIITDCYQPCAQYWVQPAQFTQHHIVTSNMIDWASCNQGHITENRWISWYWNCAHDYKRNPKFGGRDSACWFLKMTYGPMYCAFEGQVMIKDFFPRSWHKAHHVHTEHIWQRCCYKDCCLGQSMFSAQTLDWGIYRKWVGIKHVIVMFSNREMASLRRDGMHDEDQKDRCPMGIDGLTPAKTSAPASTVNQGEDILMAMKFQDCYKKFQVHWVPMEWFSLPNHYFPWYGMNPWNEELGLASHRDYIRETCPIAEEWAPIILATSWQITCEYFWVDIAQSELYPVQWQDCFYHSIAAQEVFQVPEFVFECGESLKKSGHTRWCNFDINQLMSYRFMVNMCVGADIYVQNKPYDAIECILAQAVDWWQAYLCVANTHLGAGILEALCTADCPQELDCPKAQYRCNMYIAKILIMWISNTNNNWGIYIARQCFECCNPICRAAQYEYEPSRCFFTKNWHATFNQATVCYEYHHFRVETKHMSKFIAKFTKPRDIEDVCQYQRMGRCVGRSFYYTYRHTCCYGWSWWWEMPYPCYFIHCGAAHRWCRNCVFTPLQRPGHNCDANHDCPPIYRIIIDCRQCIGLHICCKCNHHRSHPSWKNFDHHADGRHELDYRFNRDNNSGGGKCHTCHLMRIGMDGQLCHWLMSGDGYIPEFFYHNVHDDIVQHVRQVLLGIQHPYHDFEEFDNDSASVQENTTIEDSQWDLFFCYATCPDFYFPPRMRRHLGVIYSERWFDMRGINSQMEEKPQMTNHISIMMIYLEYMGICNIIINHRIKELGSQHNFTEMVRGIHTWFRTYAEPNAQWDRVWSCCYEGANNWWSMMSGATMRVKARIHHFVKIPSCIGLVHHKTTRPYISKGKVNTPYCPVDQCISAYFTGHASQQTGCWRDWHGYEVWMFSLTNYSWEGKKAWVHVTLMWDENQAWTQEQEPMPVMLMTGMFDIGRFYFVCQSRHVSNMQWLEALYFKDRCYNTMAHWEKFQSNIKMEELIGPNYVDCQYFSYIVPHSRATIENPKMMVKPQTCFPAVQYRFADIAFIMHGFHFGWAWPGINNTFPNTMYIYWRDFCIIFERWWKDGIPNGPERKPMGYDRKVQEADVMTVLQRQIKILYPIKSYYHMEFECKAYGQAHFFYGKDVAESKNNCMHDQEFHRYAGVWSTCYCFLFNFDSKRWDIIGSWEKCWITTWVTNNWKDKCQVSQHLTKHAIKSWAVDHDQFHQYTRINTNLGLHLEQYHCRCLIHCVNAQGFPYIRNSMFCCKNSILWEVLVLYEYLWKYWTSWDKLIDNDCINIECTGRGHYVHIVSCNMWIQLYDMDHHKKVKNPDIHKDDAWDKDVSCGGGEEGGPYAYILWPKMMLFHFDVKVSQYYEGNESRFHVEGLCQDGSYKPKIFSETEQCVSNDNTFSDGKYPRFINQAPQTQKVFACRFGFGAQSRHSSTAMYDNVEIPGWWNYAGWEEMQMYSYMSRFRCQCHLVTQLHSRCAVPKDQHDTGTGVMQDKSPPLVWDFALTNEPTKYHQADLLGPICGHVEGLIKLDEKRGSPTYFSEQVKMCPMFACQSINAHWCMTDATMHKVVRFYEVQYDVYERHHHHDNEDFIIAWNLKSCTQWGRVWDKDWKEHFDGFRNSIYPHKMPYNLNEPMENYPDGNLGLIQHAIQFSFPPYKIEMDSARALWMFCVNYQPMAYTEEKLFKMTCKPWVKFGIQDYRHSWASKKKGRNMHVHNADYNSFSPRESWIYCFSDCSDHGGICGLLDMMFSAVCFFCFRACSIWCKGEIIPMVPWTPRDYCMFHPLHFEYNTEDNHYPMQNFRSLFDFCLERCSEDRNEINQDYALSMATSDVYITKFPEQSFQAMKWDLRSFWAWKQMSSGDTSMPKNTGGSIKIPWDPLGVPTKSWSIHFLRDAGRGWAEAALKMEVWYPQGNNQQDWMHSCHDLPYACMWFQWEANWCKGKQPIRVCSSAGFRIRNFMTEKKNHPWHDQPCDGWTMDNIATRKNLLLGRSACAMENQKCVLVMGYNSIMGETSCGRQWFMWTDVKYAKAVIMVQAQEWIGYASMLCFKVIMKPKNAPRGDEIIPLKVNWDPRPESCGAQQMPFILWYQANIKNTAHDNWCNYFVWYTFPEHCCEDTRRGYDYANYRTTKMATKTGHGPDRYWDLMATCTGEKYNEHRKLENSHQKKTGEYGYVDEGMFSVREWQNEDLMEDSFNMEDIGVEGFQVPFTLHVVNSPKTGHQQPRDIDYCSISDPHVTHYCIVRMFSAASYVIMVWFIKHVWQAMPQALQNAPEVPYTWKYSMIGLHILKSSDNKFFCTMTTHDGRAGPVYNMTYYVKVPVNCLIWPMVMTGTVPFYQGWCSPEHKEPWIHGTASEEQDHEWHDALNTCFYQHGSCNAAVYSMGLHMGGFAIYHQFWSWQCDWDFGQQFKCDIHFIRNPPIHMYSHVDDEAVVPTRACKKYLYTHRKGVEINPYWPTTTWWYRAEDSRDFVLAQFDNMHEDHPEEMVTRLSQGMAVSWPYYEQPKTWWFRMFPPWRNSACWCQVICGCCMRLWVTIWFPESKMYCCPWHEPSQITPRNTCKPRYFAKMEITKNLTRFSMQAFVQIPNTYISAEFVPHKEQCRYGWCKYVPANSWLQQKPDYLHTSAINMMKLFNTCRTNPVNVTLPFLFWPKFNDTDVKLCPGDNCGQFHLRQWWTYPWFDSRYNQPTFLEPGDYIRASVVFQDAGEVTHVHFSDSWGKVPVMDGTAPGAATPDSMELRWCECETCIRDLFWWGPFQGYCLSTKYLGWISLQSYGPCWTFNYLERSDEQFGSAVVKMRMVDNWYIVYPGCIWCMGVHACFDQKPDNLNDNVSSSKCSVPVFMPEHQRHAMAHLFEPMMVRYLCHVSNENFETVYYSPYLQPESVTCSFHDRGMIAGVKSGRVNVHPTVKREFEKYQEHHTSSPTSEYCWPIEEPFEKPEWNTAHLREWQSANHKYETKCEPKRTWWRKSYYQDRYCYTRFLFTCIGNMNAEVFDLHCRWGMASSPNYTWMFAIPVYWTREMPFQRGHHHIRPLWSFICFHCSELFFLCKVMDEIQNKAFVDWNCCWKLYDYWGCQEPLNQETDNGGQKYYPAFRIFQKLAQATIFDQEPSVDPWYACQARMTWRTDIFSMVVPVKNNWDPTVIWYPPKRRPPHSPKTHFWFPKTMAAKVSWAQDRRRPHFNLLMGVKPLPIHDYEMDINEEFLHTTYGFWRHWEDQWQHLHTGSNKNSQRLHLFNCHSTERCHDTEKEDRDPEHASARCNHGLAHLLCKKTEEFKNPWRQHNMGQRWHQNSNLRDNQEVVKTTDVIASGYQDMYQKQHSDAHSIAGWQGELFMHMVFTMMDPPKQMRLWHTDTPFYWFWWYQKSDLAYQEWDASCFKCFTYIWGFLGHSDGVGYSCHSLWEDNICWFAKQAWAPQNGKAFDRYLCENVGPLGFEIGGKELNHNCKVKCIEDMCMVTEPINQLCPPLECHNKTGARHYKFRHPFDGPIETDGMYCCTIEKKIHHCYVG", "PEGKLHYNFSYFHEEPQVQLYCSGPDIRKVCWFYWIAGPNSKDHCTHEHLPGCCGKFAWFFRTTSREQNAKAASTMPPAMLRPEGSRQNMKMPQMFKNTRYFDYDHFQCVGKQGKNDLSIPVWKFEIDVDFPLLASANKYPWPRHTWNCYIVWHHLTNKIVEHYHLWVPLLCNTQNKPMESSLNKWPEIQLQYDCRDHVCYQYRCIRKKMIKRIFWHCWHPQHDFFATHWGAYILYGFPKCEAHNCFQQVQPAQPNQKKLDRCVREHIIFEPCDDCGCIRTHNCGRRGMPFADAGHEYTSEKVEWTYCYCRNWTMVATFYYDAICPVFNLSCEAFHPTVEMSAHGKEMIGAVMHCVMISTRQSTCMVFEFICWNDCKAHHQIFTGGFYVYRTTYRIERWWEMTYMVENQYNTWLAKAGFACYSWHEGWAHPRGGISSFPKFQPVFWDLWHAYNLKCNNPKDRDSESVKCPWSNCDMRLPKQQCVNHCTRDMRTHYSWTIPQAVYMTQWKYGAHIWHSPDKYFEEETVMYMVGNIANITRDEHCFWVLVAHASFLPRFIPFPFREDIGYHKRRYITVHDAGWQVLWNQDERITPWRFMEQHKDRIKTGSRTIRFPQSSPLTAEEQTVRMHCINCEAMMRALEQKWQCGFLSHFSYEMSIWMVKAYLPKNISFQPNWNEINGGCGAFRNLKYHTYHVQRMNKCPDGACNCHQSDDKYESMVMMWNLTYHASVTYGMVYQKQIKTGYMFMPCAFQSRPWLSDWQCDYLLPFFRYSVVPMDGHCTIPEPDPTLHFHHQTTFRGDAGAARAYIYRFHLSAKFIQIRDTVVWICEESLQGFPGYPFERAGDVDCRPHCDLFCQLCCAHAWPCDIIHLRMPVLFINFNNFMLDWYEWPEQQSPRAVQSRKTKARAMDADARNPKQRLIITDCYQPCAIYWHVLCPSNPKTIDWASCNGGKQRMALKIWWNCAHAYFSVATLKWRNPKFAGPGDSWAKSKMTYGTGTYFCMVMTINWYFFPRSWERRKAHHVHTEHIWQRCCYKDCCLGQWGIKHVIVMFSNAEMASLRRDDMHDEDQKDRCPMGIHNGLTPAKTSPNCGEDILMAMKCYFKNQGKFQVHWVVWLARMEWFSLRMKHYEHINHYFPMNPWNDELHIMFSNCNMIRETCPIGERELECFWVDIAWSELYPVQWQDCFYHSIAAQEVMQVGEESLKKIGHTRWCNFDINPSMSYRFYVQNKPYDAIECILPQAVYHKWWQALTVYKGLCVAMTHTGDSYLSMHGILEKMQVNLKSPLNQDIPQELDCPGAMYINTNANWFFEPAIYIQCFECCPAWVCNNRCFFTKNWQATFNGCAYHTVKYEYHRVETTAVHMSKFIAKFTKPEDIEDCVLTYRDYGWSWWIHCGLNAVFTPNCDANHDCPPIYRIWKLHICCHCNHHRSHCSCKNFDHHADGRFVCELDYRFNRIGTNNSGGGKCHDGQLCHWLMSGDGYIEEVDVEFTVLNWYQNFVNGQFNVAHDDVQTQHHDREEFDNDSASVQEGKCTTIEDSQNLFFCYATCPDFVFPPNENINQFCMRRHLGVIYSERWEEKPQMTGHISIMMSAECNNHRIKELGSQHNFTEMVRGIHLWFETLCVFERCNHQEHGFKEADRMMSGATRVKARIHHFVKIPSCIGLVHHKFFSYNMHSWQRPYISKGGVNCPYCPVDQCIWRCYPMCVMQSNYFTGHASQQTGCWNDKHGYEVWMFSLTWHQWGDSWRMFVPYKQWGDPAMKAWVHVTLMWDENQAWTVEQEPMPVMLMTNVPRPQGRFYFVCQSRHVSNMQWLEALQFKDRNYNTMAHWEKFQSNLRHVMFWWGPNYVWQYTPMANQLVPHSNATIENPYTFCQDHMMVKPQTCFPAVQYATAHGFHFGWAWPGFPNTMYKGYWRDNCIIFERWWKDGIPNGYERKTMGYDKCADRKVSNHISGTRAEADVMTVLYRQQVLYPSWAHWPYDNSYYHMEKPPYYECKAYGQAHMWPLVAPGVERYGKDPKNNCMHDIFATWQLGVFLFSTCYCFLFNFDSKGWDIIGSYEKFWICTWATNNWKDSCQVSQHRWHYHFRGTKHAIHQYTRINTNLGLHLEQYHCRCLDHCVNAIWMHCWTINGFPYIRNIMFCCKNSILWIMEAVLVLYEYLWYWTSWNIECTGRGHYVHIVSCPNSSVWMWIQLYDMDHHNKVKNPDIHKDDAWDKPYAYPLIPKMMLFHFDVKVSQYYEGNRSLMEIKCTKNPSIFSEAVQNTVVNTFMGSEGYPDGKYKRQRFINQAPQQKVFACRTAQSRHSGTAMYDNHEIYWATSYKEGWWNYMCRYSYMSSDVFILRYTFRIQCHLVTQLVANMALHSRCAVPKDQHGTGTGVHEMIYYTRFQDWSPPSVWTFALEMGQILLGPQCGYFMPIRIVEGLIKLDEKRGYWCFLPRYFSEQVGMKCLGMMRQPNCFACQSINFERTPITADNTMHKVVRFYEVKYDMAYQLFQGSFIIATNLKECTQWGRVWDKCWKEHFIEWPIFQKYTPGACVTSYPQAHKMLYNLNEIICEKGKMENYPDGCMHAIPFKIEMDSAWMFCLNYQPMAYRSNDKEEKLMKMTCKPGLEIQDYRHWWAGNMKKKGRNMHVHNADHPREFWIWCFSDCSDLWSYMMVGWAMCGYDMPLAECRQMKTSVPWFRACSIWCKGEIIRHFEYNTEDNHAPMRRYYPAFNCRSLFDFCLERCSEDRNEINQDYKLSMATSDVYTDATKFPEQTFQAMKWDLRYKDFMAWDQMSSGLCNNSHPTSNTGVPTKSWDIHDAGKSGWAEAALFMYPQGNNQQDWMHSMWYRKHSRQWEANVTWYCWSSAGFRIRNFMTAKQNHPWLDQPCDGVATMDNITTWHKNLLLGRSACDQFCTMENQKYNYTPQNRMGETCISHHQEAMCGREWFWWTDVKYAKLYCGHFVIMHQAQDWIGYASMTHIWQDKNAPRGDQRKVNWDPNNSPKVCKLNTTELNPHLQGSRHNIKNTAHDNVLYNAYTAAFPEHCCHDTRRPGYRTTQMAWKTGHGADRYWLMATCTGENYNENINWEHQKKTGKYGMFSVRDWQNGTQGCDLMEVIGVEGMQQVNSPKTGHQQPSEIHSREFLISDPHVSDHMMFSAAYAFGVNVIMVWQARMDWQIPCPEVPYTARYYNFCYSMTGLHILKSDNRFFCTMSHHWDMTGTYEQGRAGSVYNMTYYVGGCKMFVYEVSETPMVMTPTVPFRQGWCSPEHKFQFNTANPWIHGTRSEEQDYTHCQVEWHNTTFYQHGWSCHQVMCNAACYSMGLHMGGFAIYHQFWETETANPPIDMYSHVDDEAIDPAMSKVIGRACQKYLKGWKRGERKGVEINPYWPTTTWWYRAEVPNPMHPEEMYCNVSWPYYEQPKTWWFRSGDHFFPFWHNSACWCQVICGCCMRLWVTIWFPISPMYCCPQNTPRNTCKPWYFAKHEITKNYDDTRQSMQAFVQIQENTYISAEFQCDWFLTYYGWCKYVGPEVLYRANSWLQQKPDYLHTSAINMVWHATSTCMKLFNTCRTNRHDRNNIWGAINLFWGRPHRKNETDVVLCPGDNCGQFHLRQWWNYRDETDDVISVVNQDAGEVTHVHGRVFVPHSDSWGKVPVMMGTALRFWCNGAATRWCECEMKVNGYSCIRDLFWWGPFQGYGWISPQSYGPCWTFNYLERSWEQFGSRMVDNPIVYPGCIACMGVHACNVSSSKCSVPEHWPMMVRYLCHVSWFNFETNYYPDTEPYIVEPESVTCSFHFRGMIMWYKGMKYKSGRVSVHPTRIREEKYQEHHTSSPTSEYCWYIEEPFEIPLWNTQEPSIKHLREWQPTNHPKRTWGKLAFLFPCISNMNAEVFDLHCRWLEYMASSPNYTWKFYWTREMPFQRGHHTCICHCYSEIRPEFFITFLCQWVRELCKPMDEIQKAFVDWNCDKKEDSSCLETKLYDYWGVVNVHQWLNQETDNGGQKYYPAFTPDPIFQKVANTMWGAQAQRYMDIFDQEPSKDPWYAVEVNKARRDIFVMVKQNWDWTVKWSDHNPPGWYVSYKMRPPISPKTHFMAVSWAQDRRWGVKALPIHDHEMDINELFLYTTSKFWFTVERLHTGSNKNSYYLDSMNVYRLHPCYRDNCHSTERCHDTEKEDRQPSMDGDPEHASARCLHGLEMAGWTTNCLLCKNTEELKNPWRLHNMGQRHSHEVYFWEQNSNLRDNQEVVKTTDCTFTQKQHSHMVFIMTTPSATNMMDPPKQMNANHSDKCGKEMTDTPEEPYWFWWYQKEWDATACMHLCFKNFNCYIWGFLGNYIIMVVDIVGYSQRGVLFKMWFAYQNGDRYLCGFEGGGKELNHICKVSIYTVCMQCNTSLGLRHHHTPPLECHNKTGARHYKFRHPFDGTIETDGMYSCKIETKKHHCYSHFQG")