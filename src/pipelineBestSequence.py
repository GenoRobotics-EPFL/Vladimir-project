import sys
from dataRetriever import *
from dataCleaner import *
from consensusMedaka import getConsensus
from identification import identify
from qualityConsensus import *
from visualise import *

# todo: implement the increase in value of x if consensus doesn't work or if consensus quality is too low

# used for preprocessing
referenceReadForOrientation = None
allReads = []

minCoverageDepth = 40
x = int(1.5 * minCoverageDepth)

# at what percent in increase in consensus quality do we stop the pipeline ?
thresholdEarlyStopping = 5

outputDir = "./outputPipelineBest/"
os.system(f"rm {outputDir}*")  # clean dir to make sure we start fresh

# this is were most of the results are posted every iteration
outputFile = open(f"{outputDir}results.txt", "w+")

# this is used to plot the quality of the reads in the sample over time
allAvgQualityReads = []
# this is used to plot the quality of the consensus over time
allQualityConsensus = []


def qualityReadForSorting(read):
    """
    This method influcences the sorting of all the reads received from the start of the sequencing.
    The best reads will be at the end corresponding to the choice of features for the meaning of "best read".

    the quality of the read can be influenced by:
    - its length
    - its mean quality score (the higher the better)
    - its mode quality score (the higher the better)
    """

    return len(read.sequence) * getQualityRead(read)


def getReadsFirstIteration(geneName):

    # retrieve new reads
    reads = waitForReadsFromFile("fastqpass/iteration0.fastq")

    global referenceReadForOrientation
    cleanReads, referenceReadForOrientation = getCleanReads(
        reads, geneName)

    global allReads
    allReads = cleanReads

    allReads.sort(key=qualityReadForSorting)


def getReadsLaterIterations(geneName, iterationNum):

    # retrieve new reads
    reads = waitForReadsFromFile(f"fastqpass/iteration{iterationNum}.fastq")

    cleanReads, _ = getCleanReads(
        reads, geneName, referenceReadForOrientation)

    # append new reads at the front as they are more likely to have a low quality
    global allReads
    allReads += cleanReads

    allReads.sort(key=qualityReadForSorting)


def getBestXReads():
    # bestXReads = allCleanReadsSorted[:]

    # tenHighest = int(len(allCleanReadsSorted) * 0.05)
    # bestXReads = bestXReads[:-tenHighest]

    bestXReads = allReads[-x:]

    return bestXReads


def createNewConsensus(bestXReads):

    # write those reads to file so consensus can act on those reads
    pathReads = f"{outputDir}temp_readsForConsensus.fasta"
    writeReadsToFile(pathReads, bestXReads)

    consensus = getConsensus(pathReads)

    return consensus


def saveQuality(bestXReads, iterationNum):

    qualityConsensus, coverages = getQualityConsensus(
        "./outputMedaka/calls_to_draft.bam")
    allQualityConsensus.append(qualityConsensus)

    saveCoverageVisualisation(outputDir,
                              qualityConsensus, coverages, x, iterationNum)

    qualityReads = sum(map(qualityReadForSorting, bestXReads)) / x
    allAvgQualityReads.append(qualityReads)

    saveQualityReadsOverTime(outputDir, allAvgQualityReads)
    saveQualityConsensusOverTime(outputDir, allQualityConsensus)
    saveQualityBothOverTime(outputDir, allAvgQualityReads, allQualityConsensus)


def getIdentification(consensus, db):

    pathCons = f"{outputDir}tempConsForIdent.fasta"
    writeConsensus(pathCons, consensus)

    # db = "rbcL" # "psbA-trnH"
    name, cov, iden = identify(pathCons, db)

    return name, cov, iden


def checkEarlyStoppingCriteria(iterationNum):

    if iterationNum <= 3:
        return False

    x1 = iterationNum - 4
    x2 = iterationNum
    y1 = allQualityConsensus[x1]
    y2 = allQualityConsensus[x2]

    maxQual = max(allQualityConsensus)
    minQual = min(allQualityConsensus)

    percentIncreaseSinceStart = ((y2 - y1) / (maxQual - minQual)) * 100

    return percentIncreaseSinceStart < thresholdEarlyStopping


def start(geneName):

    print(f"Starting best x pipeline with gene name {geneName}")

    iterationNum = 0

    while True:

        print(f" --- Starting iteration: {iterationNum} ---")

        # get the new reads
        if iterationNum == 0:
            getReadsFirstIteration(geneName)
        else:
            getReadsLaterIterations(geneName, iterationNum)

        # take the best x ones received so far
        bestXReads = getBestXReads()

        # create consensus based on those
        consensus = createNewConsensus(bestXReads)

        # save data from this iteratoin
        saveQuality(bestXReads, iterationNum)

        # get identification
        name, cov, iden = getIdentification(consensus, geneName)

        # print results to file
        outputFile.write(f"iteration: {iterationNum}\n")
        outputFile.write(f"name entry: {name}\n")
        outputFile.write(f"coverage: {cov}\n")
        outputFile.write(f"identity: {iden}\n")
        outputFile.write(f"consensus: \n{consensus}\n")
        outputFile.write("\n")
        outputFile.flush()

        # check if you can stop the pipeline
        canStop = checkEarlyStoppingCriteria(iterationNum)
        if canStop:
            print("Early Stopping criterion met. Stopping pipeline")
            break

        iterationNum += 1


def main():

    geneName = sys.argv[1]

    start(geneName)


if __name__ == "__main__":
    main()
