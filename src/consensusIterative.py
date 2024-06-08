"""
This file was an attempt to create a iterative metho for creating a consensus where all the reads are aligned on the previous consensus.
It didn't work.
"""

import sys
import os
import spoa
import matplotlib.pyplot as plt

from utils import *
from dataRetriever import *


def spoaCons(sequences):
    """
    get a conensus of multiple sequences (bases only), using the spoa library.
    the consensus will compare all sequences together to find the best match.
    """

    # run spoa
    min_cov = int(round(len(sequences) * 0.15))
    cons, _ = spoa.poa(sequences, min_coverage=min_cov,
                       genmsa=False)  # , n=0, g=-127, e=-127)

    # run spoa a second time with the previous result as first read
    # , n=0, g=-127, e=-127)
    cons, msa = spoa.poa([cons, *sequences], min_coverage=min_cov, genmsa=True)

    return cons, msa[1:]  # first line of msa is the cons we added


class Score:
    allowedBases = set(['A', 'T', 'G', 'C'])

    def __init__(self):
        self.frequencies = {'A': 0, 'T': 0, 'C': 0, 'G': 0}
        self.bestQScores = {'A': 0, 'T': 0, 'C': 0, 'G': 0}

    def totalCoverage(self):
        return sum(self.frequencies.values())

    def acceptBase(self, base, qscoreChar="0"):
        if base not in self.allowedBases:  # i.e base == '-'
            return

        qscore = ord(qscoreChar) - 33

        self.frequencies[base] += 1
        self.bestQScores[base] = max(self.bestQScores[base], qscore)

    def getBestBase(self):

        # best base is the one that appears the most often
        highestFrequencyBase = max(
            self.frequencies.items(), key=(lambda item: item[1]))[0]

        # best base is the one with best score
        highestQScoreBase = max(self.bestQScores.items(),
                                key=(lambda item: item[1]))[0]

        return highestFrequencyBase

    def __str__(self) -> str:
        return str(self.frequencies)


class TemporaryConsensus:
    """
    this class is created once with the very initial conensus of the first iteration.

    the reads of the first iteration and the subsequent ones will then be alligned with it
    """

    def __init__(self, reads, msa):

        self.numReads = len(reads)
        self.sequence = []
        self.scores = []

        baseIdxInRead = [0 for _ in range(self.numReads)]

        for baseIdx in range(len(msa[0])):

            score = Score()

            for readIdx in range(self.numReads):

                base = msa[readIdx][baseIdx]

                if base == '-':
                    continue

                score.acceptBase(
                    base, reads[readIdx].quality[baseIdxInRead[readIdx]])

                baseIdxInRead[readIdx] += 1

            self.scores.append(score)
            self.sequence.append(score.getBestBase())

        # self.trimSequenceEnds(numReads)
        # self.removeLowCoverage()

    def removeLowCoverage(self):

        thresholdFrequency = self.numReads * 0.15

        newSequence = []
        newScores = []

        for i in range(len(self.sequence)):

            bestBase = self.sequence[i]
            frequency = self.scores[i].frequencies[bestBase]

            if frequency >= thresholdFrequency:
                newSequence.append(self.sequence[i])
                newScores.append(self.scores[i])

        self.sequence = newSequence
        self.scores = newScores

    def trimSequenceEnds(self):

        thresholdCoverage = 0.05 * self.numReads

        # trim the end of sequence
        numBasesConsensus = len(self.sequence)
        startOfWalkForward = int(0.6 * numBasesConsensus)
        startOfWalkBackwards = int(0.4 * numBasesConsensus)

        for i in range(startOfWalkForward, numBasesConsensus - 2):

            if self.scores[i].totalCoverage() <= thresholdCoverage and \
               self.scores[i + 1].totalCoverage() <= thresholdCoverage and \
               self.scores[i + 2].totalCoverage() <= thresholdCoverage:

                self.scores = self.scores[:i]
                self.sequence = self.sequence[:i]

                break

        # trim the start of sequence
        for i in range(startOfWalkBackwards, 2, -1):

            if self.scores[i].totalCoverage() <= thresholdCoverage and \
               self.scores[i - 1].totalCoverage() <= thresholdCoverage and \
               self.scores[i - 2].totalCoverage() <= thresholdCoverage:

                self.scores = self.scores[i:]
                self.sequence = self.sequence[i:]

                break

    def removeExtentionOfSequence(self, msa):

        # trim from the start
        newStart = 0
        while msa[0][newStart] == "-":
            newStart += 1

        # trim from the end
        newEnd = len(msa[0]) - 1
        while msa[0][newEnd] == "-":
            newEnd -= 1

        newMsa = [row[newStart:newEnd] for row in msa]

        return newMsa

    def acceptRead(self, newRead):
        self.numReads += 1

        # the logic is that we don't want a new read to extend the current sequence too much
        # as we know that we have already found a good consensus with good enough length

        seqAsString = "".join(self.sequence)
        _, msa = spoa.poa([seqAsString, newRead.sequence], n=-1)

        # msa = self.removeExtentionOfSequence(msa)

        # writeMsaToFile("output.fata", msa)

        self.mergeAlignment(msa)

        # once in a while, trim the sequence
        # if self.numReads % 50 == 0:
        #    self.removeLowCoverage()

    def mergeAlignment(self, alignment):

        sequenceIdx = 0

        for i in range(len(alignment[0])):

            baseCons = alignment[0][i]
            baseNewRead = alignment[1][i]

            if baseCons == '-':
                continue

            elif baseNewRead == '-':
                sequenceIdx += 1
                continue

            else:
                self.scores[sequenceIdx].acceptBase(baseNewRead)
                self.sequence[sequenceIdx] = self.scores[sequenceIdx].getBestBase()

                sequenceIdx += 1

    def getCoverageSequence(self):

        depths = []

        for score in self.scores:

            depths.append(score.totalCoverage())

        return depths

    def getCons(self):
        """
        basically self.sequence contains the whole messed up alignment.
        here we clean it to find the correct consensus
        """

        thresholdCoverage = 0.15 * self.numReads

        polishedCons = ""

        for i in range(len(self.sequence)):

            bestBase = self.sequence[i]

            totalCoverage = self.scores[i].totalCoverage()

            frequencyBestBase = self.scores[i].frequencies[bestBase]
            scoreBestBase = self.scores[i].bestQScores[bestBase]

            # and frequencyBestBase >= (0.30 * totalCoverage):
            if totalCoverage >= thresholdCoverage:
                polishedCons += bestBase

        return polishedCons


def plotmsadepth(obj):

    plt.plot(obj.getCoverageSequence())
    plt.show()


def main():

    pathfastq = os.path.join(sys.argv[1])  # this must be the cleaned reads
    reads = getReadsFromFile(pathfastq)

    reads = reads[:10]

    sequences = [read.sequence for read in reads]
    consOfSpoa, msa = spoaCons(sequences)

    tempConsObj = TemporaryConsensus(reads, msa)

    originalCons = tempConsObj.getCons()
    print(f"Original consensus: \n{originalCons}")

    # # improve consensus with new reads

    # pathfastq = os.path.join(sys.argv[2]) # this must be the cleaned reads
    # newReads = getReadsFromFile(pathfastq)

    # for newRead in newReads:
    #     tempConsObj.acceptRead(newRead)

    # print(f"Consensus improved: \n{tempConsObj.getCons()}")

    # plotmsadepth(tempConsObj)


if __name__ == "__main__":
    main()
