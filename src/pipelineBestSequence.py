import sys
from dataRetriever import *
from dataCleaner import *
from consensusMedaka import getConsensus
from identification import identify

maxLength = 750
minReadQ = 10 

allCleanReads = []
referenceReadForOrientation = None

consensus = None

outputFile = open("pipelineBestOutput.txt", "w+")


def getReadsFirstIteration():

    # retrieve new reads
    reads = getReadsFromFile("fastqpass/iteration0.fastq")

    global referenceReadForOrientation
    cleanReads, referenceReadForOrientation = getCleanReads(
        reads, maxLength, minReadQ)

    global allCleanReads
    allCleanReads += cleanReads


def getReadsLaterIterations(iterationNum):

    # retrieve new reads
    reads = getReadsFromFile(f"fastqpass/iteration{iterationNum}.fastq")

    cleanReads, _ = getCleanReads(
        reads, maxLength, minReadQ, referenceReadForOrientation)

    global allCleanReads
    allCleanReads += cleanReads


def getIdentification(iterationNum):

    pathCons = "tempCons.fasta"
    writeConsensus(pathCons, consensus)

    name, cov, iden = identify(pathCons, "psbA-trnH")

    outputFile.write(f"iteration: {iterationNum}\n")
    outputFile.write(f"name: {name}\n")
    outputFile.write(f"coverage: {cov}\n")
    outputFile.write(f"identity: {iden}\n")
    outputFile.write(f"consensus: \n{consensus}\n")
    outputFile.write("\n")

    

def updateNewConsensus():


    bestReads = sorted(allCleanReads, key=lambda read: getQualityRead(read))

    # remove top 10% of best quality reads
    # todo: instead of removing the 10 percent, remove the ones that have reads with a quality too high
    bestReads = bestReads[:-int(len(bestReads) * 0.10)]
    # only keep 100 reads at a maximuml
    bestReads = bestReads[-100:]

    # write those reads to file so consensus can act on those
    pathReads = "temp.fasta"
    writeReadsToFile(pathReads, bestReads)

    global consensus
    consensus = getConsensus(pathReads)


def start():

    iterationNum = 0

    while True:

        print(f"Starting iteration: {iterationNum}")

        if iterationNum == 0:
            getReadsFirstIteration()

        else:
            getReadsLaterIterations(iterationNum)

        updateNewConsensus()

        getIdentification(iterationNum)

        iterationNum += 1
        outputFile.flush()

    outputFile.close()


def main():
    start()

if __name__ == "__main__":
    main()
