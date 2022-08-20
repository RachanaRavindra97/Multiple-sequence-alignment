import time
from sys import getsizeof
import math as math
import numpy as np

MATCH = 5
MISMATCH = -4
INDEL = -8


def create_scoring_matrix():
    i = 0
    j = 0
    k = 0

    scoringArray = np.zeros((21, 21, 21))

    for i in range(0, 21):
        for j in range(0, 21):
            for k in range(0, 21):
                score = 0
                if (i == 20) and (j == 20) and (k == 20):
                    score = 0
                elif (i == 20 and j == 20 and k != 20) or (i == 20 and j != 20 and k == 20) or (i != 20 and j == 20 and k == 20):
                    score = INDEL * 2
                elif (i == 20 and j != 20 and k != 20):
                    if (j == k):
                        score = (INDEL * 2) + MATCH
                    else:
                        score = (INDEL * 2) + MISMATCH
                elif (i != 20 and j == 20 and k != 20):
                    if (i == k):
                        score = (INDEL * 2) + MATCH
                    else:
                        score = (INDEL * 2) + MISMATCH
                elif (i != 20 and j != 20 and k== 20):
                    if (i == j):
                        score = (INDEL * 2) + MATCH
                    else:
                        score = (INDEL * 2) + MISMATCH
                elif ( i < 20 and j < 20 and k < 20):
                    if (i == j and i == k):
                        score = MATCH * 3
                    elif ( i == j and i != k) or (i != j and i == k) or (j == k and j != i) or (j != k and j == i):
                        score = MATCH + (2 * MISMATCH)
                    elif (i != j and j != k and k != i):
                        score = MISMATCH * 3

                scoringArray[i, j, k] = score

    return scoringArray


def alignment_seq(iStart, jStart, kStart, lengthSeq1, lengthSeq2, lengthSeq3):

    # Score Forward

    lengthi = lengthSeq1 - iStart
    lengthj = lengthSeq2 - jStart
    lengthk = lengthSeq3 - kStart

    print (lengthi)
    if (lengthi <= 1):
        return

    mid = math.floor((lengthSeq1 + iStart) / 2)

    score1 = np.zeros((2, int(lengthj + 1), int(lengthk + 1)))
    for i in range(0, 2):
        for j in range(0, int(lengthj + 1)):
            for k in range(0, int(lengthk + 1)):
                score1[i, j, k] = -100000
    score1[0, 0, 0] = 0

    for k in range(1, lengthk + 1):
        score1[0, 0, k] = scoredArray[4, 4, seqLen3Arr[kStart + k]] * k
    for j in range(1, lengthj + 1):
        score1[0, j, 0] = scoredArray[4, seqLen2Arr[jStart + j], 4]

    for j in range(0, lengthj + 1):
        for k in range(0, lengthk + 1):
            score_xy_plane = 0
            subscore = np.array([-10000, -10000, -10000])

            subscore[0] = score1[0, j - 1, k - 1] + scoredArray[4, seqLen2Arr[jStart + j], seqLen3Arr[kStart + k]]
            subscore[1] = score1[0, j - 1, k] + scoredArray[4, seqLen2Arr[jStart + j], 4]
            subscore[2] = score1[0, j, k - 1] + scoredArray[4, 4, seqLen3Arr[kStart + k]]

            score_xy_plane = subscore[0]

            for i in range(1, 3):
                if subscore[i] > score_xy_plane:
                    score_xy_plane = subscore[i]

            score1[0, j, k] = score_xy_plane

    lengthy = 0

    for i in range(1, (mid - iStart + 1)):
        lengthy = (lengthy + 1) % 2
        score1[lengthy, 0, 0] = score1[(lengthy + 1) % 2, 0, 0] + scoredArray[seqLen1Arr[iStart + i], 4, 4]
        for x in range(1, lengthj + 1):
            score_xy_plane = 0
            subscore2 = np.array([-10000, -10000, -10000])
            subscore2[0] = score1[(lengthy + 1)%2, x - 1, 0] + scoredArray[seqLen1Arr[iStart + i], seqLen2Arr[jStart + x], 4]
            subscore2[1] = score1[(lengthy + 1)%2, x, 0] + scoredArray[seqLen1Arr[iStart + i], 4, 4]
            subscore2[2] = score1[lengthy, x - 1, 0] + scoredArray[4, seqLen2Arr[jStart + x], 4]
            score_xy_plane = subscore2[0]

            for j in range(1, 3):
                if subscore2[j] > score_xy_plane:
                    score_xy_plane = subscore2[j]

            score1[lengthy, x, 0] = score_xy_plane

        for z in range(1, lengthk + 1):

            subscore3 = np.array([-10000, -10000, -10000])
            score_xy_plane = 0
            subscore3[0] = score1[((lengthy + 1) % 2), 0, z - 1] + scoredArray[seqLen1Arr[iStart + i], 4, seqLen3Arr[kStart + z]]
            subscore3[1] = score1[((lengthy + 1) % 2), 0, z] + scoredArray[seqLen1Arr[iStart + i], 4, 4]
            subscore3[2] = score1[lengthy, 0, z - 1] + scoredArray[4, 4, seqLen3Arr[kStart + z]]

            score_xy_plane = subscore3[0]

            for j in range(1, 3):
                if subscore3[j] > score_xy_plane:
                    score_xy_plane = subscore3[j]

            score1[lengthy, 0, z] = score_xy_plane

        for j in range(1, lengthj):
            for k in range(1, lengthk):
                subscore4 = np.array([-10000, -10000, -10000, -10000, -10000, -10000, -10000])
                score_xy_plane = 0

                subscore4[0] = score1[(lengthy + 1)%2, j - 1, k - 1] + scoredArray[seqLen1Arr[iStart + i], seqLen2Arr[jStart + j], seqLen3Arr[kStart + k]]
                subscore4[1] = score1[(lengthy + 1)%2, j - 1, k] + scoredArray[seqLen1Arr[iStart + i], seqLen2Arr[jStart + j], 4]
                subscore4[2] = score1[(lengthy + 1)%2, j, k -1] + scoredArray[seqLen1Arr[iStart + i], 4, seqLen3Arr[kStart + k]]
                subscore4[3] = score1[lengthy, j -1, k - 1] + scoredArray[4, seqLen2Arr[jStart +j], seqLen3Arr[kStart + k]]
                subscore4[4] = score1[(lengthy + 1)%2, j, k] + scoredArray[seqLen1Arr[iStart + i], 4, 4]
                subscore4[5] = score1[lengthy, j - 1, k] + scoredArray[4, seqLen2Arr[jStart + j], 4]
                subscore4[6] = score1[lengthy, j, k - 1] + scoredArray[4, 4, seqLen3Arr[kStart + k]]

                score_xy_plane = subscore4[0]

                for w in range(1, 7):
                    if subscore4[w] > score_xy_plane:
                        score_xy_plane = subscore4[w]

                score1[lengthy, j, k] = score_xy_plane

    forwardScore = lengthy

    #Score Backward

    score2 = np.zeros((2, int(lengthj + 1), int(lengthk + 1)))

    score2[0, lengthj, lengthk] = 0
    for k in range(lengthk, 0, -1):
        score2[0, lengthj, k] = scoredArray[4, 4, seqLen3Arr[kStart + k]] * (lengthk - k)

    for j in range(lengthj, 0,-1):
        score2[0, j, lengthk] = scoredArray[4, seqLen2Arr[jStart + j], 4] * (lengthj - j)

    for j in range(lengthj - 1, 0, -1):
        for k in range(lengthk - 1, 0, -1):
            score_xy_plane = 0
            subscore5 = [-10000, -10000, -10000]

            subscore5[0] = score2[0, j + 1, k + 1] + scoredArray[4, seqLen2Arr[jStart + j + 1], seqLen3Arr[kStart + k + 1]]
            subscore5[1] = score2[0, j + 1, k] + scoredArray[4, seqLen2Arr[jStart + j + 1], 4]
            subscore5[2] = score2[0, j, k + 1] + scoredArray[4, 4, seqLen3Arr[kStart + k + 1]]

            score_xy_plane = subscore5[0]

            for i in range(1, 3):
                if subscore5[i] > score_xy_plane:
                    score_xy_plane = subscore5[i]

            score2[0, j, k] = score_xy_plane

    lengthy = 0

    for i in range(lengthi - 1, mid - iStart - 1, -1):

        lengthy = (lengthy + 1) % 2

        score2[lengthy, lengthj, lengthk] = score2[(lengthy + 1) % 2, lengthj, lengthk] + scoredArray[seqLen1Arr[iStart + i + 1], 4, 4]
        for x in range(lengthj - 1, -1, -1):
            subscore6 = np.array([-10000, -10000, -10000])
            score_xy_plane = 0

            subscore6[0] = score2[(lengthy + 1) % 2, x + 1, lengthk] + scoredArray[seqLen1Arr[iStart + i + 1], seqLen2Arr[jStart + x + 1], 4]
            subscore6[1] = score2[(lengthy + 1) % 2, x, lengthk] + scoredArray[seqLen1Arr[iStart + i + 1], 4, 4]
            subscore6[2] = score2[lengthy, x + 1, lengthk] + scoredArray[4, seqLen2Arr[jStart + x + 1], 4]

            score_xy_plane = subscore6[0]

            for i in range(1, 3):
                if subscore6[i] > score_xy_plane:
                    score_xy_plane = subscore6[i]


            score2[lengthy, x, lengthk] = score_xy_plane

        for x in range(lengthk - 1, -1, -1):
            subscore7 = np.array([-10000, -10000, -10000])
            score_xy_plane = 0

            subscore7[0] = score2[(lengthy + 1) % 2, lengthj, x + 1] + scoredArray[seqLen1Arr[iStart + i + 1], 4, seqLen3Arr[kStart + x + 1]]
            subscore7[1] = score2[(lengthy + 1) % 2, lengthj, x] + scoredArray[seqLen1Arr[iStart + i + 1], 4, 4]
            subscore7[2] = score2[lengthy, lengthj, x + 1] + scoredArray[4, 4, seqLen3Arr[kStart + x + 1]]

            score_xy_plane = subscore7[0]

            for i in range(1, 3):
                if subscore7[i] > score_xy_plane:
                    score_xy_plane = subscore7[i]

            score2[lengthy, lengthj, x] = score_xy_plane

        for j in range(lengthj - 1, -1, -1):
            for k in range(lengthk - 1, -1, -1):

                subscore8 = np.array([-10000, -10000, -10000, -10000, -10000, -10000, -10000])

                score_xy_plane = 0

                subscore8[0] = score2[(lengthy + 1) % 2, j + 1, k + 1] + scoredArray[seqLen1Arr[iStart + i + 1], seqLen2Arr[jStart + j + 1], seqLen3Arr[kStart + k + 1]]
                subscore8[1] = score2[(lengthy + 1) % 2, j + 1, k] + scoredArray[seqLen1Arr[iStart + i + 1], seqLen2Arr[jStart + j + 1], 4]
                subscore8[2] = score2[(lengthy + 1) % 2, j, k + 1] + scoredArray[seqLen1Arr[iStart + i + 1], 4, seqLen3Arr[kStart + k + 1]]
                subscore8[3] = score2[lengthy, j + 1, k + 1] + scoredArray[4, seqLen2Arr[jStart + j + 1], seqLen3Arr[kStart + k + 1]]
                subscore8[4] = score2[(lengthy + 1) % 2, j + 1, k] + scoredArray[seqLen1Arr[iStart + i + 1], 4, 4]
                subscore8[5] = score2[lengthy, j + 1, k] + scoredArray[4, seqLen2Arr[jStart + j + 1], 4]
                subscore8[6] = score2[lengthy, j, k + 1] + scoredArray[4, 4, seqLen3Arr[kStart + k + 1]]

                score_xy_plane = subscore8[0]

                for i in range(1, 7):
                    if subscore8[i] > score_xy_plane:
                        score_xy_plane = subscore8[i]

                score2[lengthy, j, k] = score_xy_plane

    backwardScore = lengthy

    maxScore = score1[forwardScore, 0 , 0] + score2[backwardScore, 0, 0]

    print ("mid: " + str(mid))
    midX = mid
    midY = jStart
    midZ = kStart

    for j in range(0, lengthj):
        for k in range(0, lengthk):
            if (score1[forwardScore, j, k] + score2[backwardScore, j, k]) > maxScore:
                maxScore = score1[forwardScore, j, k] + score2[backwardScore, j, k]
                midY = jStart + j
                midZ = kStart + k

    global maxValue

    if maxValue:
        maxValue = False
        print ("Final Alignment Score: " + str(maxScore))

    alignment_seq(iStart, jStart, kStart, midX, midY, midZ)
    alignment_seq(midX, midY, midZ, lengthSeq1, lengthSeq2, lengthSeq3)

    return

def protein_seq_mapping(character):

    if str(character) == "A":
        return 0
    elif str(character) == "C":
        return 1
    elif str(character) == "D":
        return 2
    elif str(character) == "E":
        return 3
    elif str(character) == "F":
        return 4
    elif str(character) == "G":
        return 5
    elif str(character) == "H":
        return 6
    elif str(character) == "I":
        return 7
    elif str(character) == "K":
        return 8
    elif str(character) == "L":
        return 9
    elif str(character) == "M":
        return 10
    elif str(character) == "N":
        return 11
    elif str(character) == "P":
        return 12
    elif str(character) == "Q":
        return 13
    elif str(character) == "R":
        return 14
    elif str(character) == "S":
        return 15
    elif str(character) == "T":
        return 16
    elif str(character) == "V":
        return 17
    elif str(character) == "W":
        return 18
    elif str(character) == "Y":
        return 19
    else:
        return 20


def main():
    global scoredArray
    scoredArray = create_scoring_matrix()

    global maxValue
    maxValue = True

    seq1 = open("seq1.txt", "r")
    seq2 = open("seq2.txt", "r")
    seq3 = open("seq3.txt", "r")

    seqLen1 = []
    seqLen2 = []
    seqLen3 = []

    for i in seq1.readlines():
        for letter in i:
            seqLen1.append(protein_seq_mapping(letter))

    global seqLen1Arr
    seqLen1Arr = np.array(seqLen1)

    for i in seq2.readlines():
        for letter in i:
            seqLen2.append(protein_seq_mapping(letter))

    global seqLen2Arr
    seqLen2Arr = np.array(seqLen2)

    for i in seq3.readlines():
        for letter in i:
            seqLen3.append(protein_seq_mapping(letter))

    global seqLen3Arr
    seqLen3Arr = np.array(seqLen3)

    alignment_seq(0, 0, 0, len(seqLen1Arr) - 1, len(seqLen2Arr) - 1, len(seqLen3Arr) - 1)

main()