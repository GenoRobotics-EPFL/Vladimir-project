"""
the idea of this pipeline is to generate a consensus with medaka on the first iteration.
On the following iterations, we use medaka to create a consensus where we insert the previous consensus as a read.

As the result of medaka is a consensus with quality, it was though the quality of the consensus should always improve

problem: medaka doesn't increase the quality of the consensus over time. On the contrary, it decreases over time, which isn't what we want
"""

import sys
from dataRetriever import *
from dataCleaner import *
from consensusMedaka import *
from identification import identify


maxNumBasesInRead = 750
minReadQscore = 10

referenceReadForOrientation = None

consensus = None

outputFile = open("pipelineUpgradeOutput.txt", "w+")


def getReadsFirstIteration(geneName):

    # retrieve new reads
    reads = getReadsFromFile("fastqpass/iteration0.fastq")

    global referenceReadForOrientation
    cleanReads, referenceReadForOrientation = getCleanReads(
        reads, geneName, minReadQscore)

    return cleanReads


def getReadsLaterIterations(geneName, iterationNum):

    # retrieve new reads
    reads = getReadsFromFile(f"fastqpass/iteration{iterationNum}.fastq")

    cleanReads, _ = getCleanReads(
        reads, geneName, minReadQscore, referenceReadForOrientation)

    return cleanReads


def getIdentification(db):

    pathCons = "tempCons.fasta"
    writeConsensus(pathCons, consensus.sequence)

    # db = "rbcL" # "psbA-trnH"
    name, cov, iden = identify(pathCons, db)

    return name, cov, iden


def createInitialConsensus(reads):

    # write those reads to file so consensus can act on those
    pathReads = "temp.fasta"
    writeReadsToFile(pathReads, reads)

    global consensus
    consensus = getConsensusAndQuality(pathReads)


def updateNewConsensus(reads):
    global consensus

    # write those reads to file so consensus can act on those
    pathReads = "temp.fasta"
    writeReadsToFile(pathReads, [consensus] + reads)

    consensus = getConsensusAndQuality(pathReads)


def start(db):

    iterationNum = 0

    while True:

        print(f" --- Starting iteration: {iterationNum} ---")

        # get reads
        if iterationNum == 0:

            reads = getReadsFirstIteration(db)

            createInitialConsensus(reads)

        else:
            reads = getReadsLaterIterations(db, iterationNum)

            updateNewConsensus(reads)

        print(f"Quality Consensus: {getQualityRead(consensus)}")

        # get identification
        name, cov, iden = getIdentification(db)

        # print results
        outputFile.write(f"iteration: {iterationNum}\n")
        outputFile.write(f"name: {name}\n")
        outputFile.write(f"coverage: {cov}\n")
        outputFile.write(f"identity: {iden}\n")
        outputFile.write(f"consensus: \n{consensus}\n")
        outputFile.write("\n")
        outputFile.flush()

        iterationNum += 1

    outputFile.close()


def main():

    db = sys.argv[1]

    start(db)


if __name__ == "__main__":
    main()
