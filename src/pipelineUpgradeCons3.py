import sys
from dataRetriever import *
from dataCleaner import *
from consensusMedaka import *
from consensusMinimap import *
from identification import identify

"""
the idea of this pipeline is to generate a consensus with medaka on the first iteration.
On the following iterations, we use minimap to map the new reads on the previous consensus, and we use racon to generate the final consensus of the iteration.
"""

minReadQscore = 10
referenceReadForOrientation = None
consensus = None

outputFile = open("pipelineUpgrade3Output.txt", "w+")


def getReadsFirstIteration(gene):

    # retrieve new reads
    reads = getReadsFromFile("fastqpass/iteration0.fastq")

    global referenceReadForOrientation
    cleanReads, referenceReadForOrientation = getCleanReads(
        reads, gene, minReadQscore)

    return cleanReads


def getReadsLaterIterations(gene, iterationNum):

    # retrieve new reads
    reads = getReadsFromFile(f"fastqpass/iteration{iterationNum}.fastq")

    cleanReads, _ = getCleanReads(
        reads, gene, minReadQscore, referenceReadForOrientation)

    return cleanReads


def getIdentification(geneName):

    pathCons = "tempCons.fasta"
    writeConsensus(pathCons, consensus)

    # db = "rbcL" # "psbA-trnH"
    name, cov, iden = identify(pathCons, geneName)

    return name, cov, iden


def createInitialConsensus(reads):

    # write those reads to file so consensus can act on those
    pathReads = "temp-newReadsIteration.fasta"
    writeReadsToFile(pathReads, reads)

    global consensus
    consensus = getConsensus(pathReads)


def updateNewConsensus(reads):
    global consensus

    # write those reads to file so consensus can act on those
    pathReads = "temp-newReadsIteration.fastq"
    writeReadsToFile(pathReads, reads)

    pathToDraft = "tempCons.fasta"
    writeConsensus(pathToDraft, consensus)

    consensus = getConsensusMinimap(pathReads, pathToDraft)


def start(geneName):

    iterationNum = 0

    while True:

        print(f" --- Starting iteration: {iterationNum} ---")

        # get reads
        if iterationNum == 0:

            reads = getReadsFirstIteration(geneName)

            createInitialConsensus(reads)

        else:
            reads = getReadsLaterIterations(geneName, iterationNum)

            updateNewConsensus(reads)

        # get identification
        name, cov, iden = getIdentification(geneName)

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

    geneName = sys.argv[1]

    start(geneName)


if __name__ == "__main__":
    main()
