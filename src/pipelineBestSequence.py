"""
This script is the best-x pipeline. See the README of the repo to see how to use it.
"""


import sys
from dataRetriever import *
from dataCleaner import *
from consensusMedaka import *
from identification import identify
from qualityConsensus import *
from visualise import *
import os

minCoverageDepth = 40
x = int(1.5 * minCoverageDepth)

# at what percent in increase in consensus quality do we stop the pipeline ?
thresholdEarlyStopping = 5

outputDir = "./outputPipelineBest/"
if not os.path.exists(outputDir):
    os.makedirs(outputDir)
else:
    os.system(f"rm {outputDir}*")


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


def getReadsIteration(geneName, iterationNum, referenceReadForOrientation):

    # retrieve new reads
    reads = waitForReadsFromFile(f"fastqpass/iteration{iterationNum}.fastq")

    cleanReads, referenceReadForOrientationNew = getCleanReads(
        reads, geneName, referenceReadForOrientation)

    return cleanReads, referenceReadForOrientationNew


def createNewConsensus(sampleReads):

    # write those reads to file so consensus can act on those reads
    pathReads = f"{outputDir}temp_readsForConsensus.fasta"
    writeReadsToFile(pathReads, sampleReads)

    try:
        consensus = getConsensus(pathReads)
    except NoConsensusEx:
        return "", False

    return consensus, True


def saveQuality(iterationNum, allQualityConsensus, allAvgQualityReads, sampleReads):

    qualityConsensus, coverages = getQualityConsensus(
        "./outputMedaka/calls_to_draft.bam")
    allQualityConsensus.append(qualityConsensus)

    saveCoverageVisualisation( outputDir, qualityConsensus, coverages, x, iterationNum)

    avgQualityReads = sum(
        map(qualityReadForSorting, sampleReads))
    allAvgQualityReads.append(avgQualityReads)

    saveQualityReadsOverTime(outputDir, allAvgQualityReads)
    saveQualityConsensusOverTime(outputDir, allQualityConsensus)
    saveQualityBothOverTime(outputDir, allAvgQualityReads, allQualityConsensus)

    return coverages


def getIdentification(consensus, db):

    pathCons = f"{outputDir}tempConsForIdent.fasta"
    writeConsensus(pathCons, consensus)

    name, cov, iden = identify(pathCons, db)

    return name, cov, iden


def checkEarlyStoppingCriteria(iterationNum, allQualityConsensus, coverages):

    if iterationNum <= 6:
        return False

    x1 = iterationNum - 7
    x2 = iterationNum - 6
    meanBefore = (
        allQualityConsensus[x1] + allQualityConsensus[x2]) / 2

    x1 = iterationNum
    x2 = iterationNum - 1
    meanAfter = (
        allQualityConsensus[x1] + allQualityConsensus[x2]) / 2

    maxQual = max(allQualityConsensus)
    minQual = min(allQualityConsensus)

    percentIncreaseSinceStart = (
        (meanAfter - meanBefore) / (maxQual - minQual)) * 100

    firstCondition = percentIncreaseSinceStart < thresholdEarlyStopping

    secondCondition = np.median(coverages) >= minCoverageDepth

    return firstCondition and secondCondition


def start(geneName):
    global x

    referenceReadForOrientation = None
    sampleReads = []

    # this is were most of the results are posted every iteration
    outputFile = open(f"{outputDir}result.txt", "w+")

    # this is used to plot the quality of the reads in the sample over time
    allAvgQualityReads = []
    # this is used to plot the quality of the consensus over time
    allQualityConsensus = []

    timeStart = []
    timeAfterGettingReads = []
    timeAfterConsensus = []
    timeAfterQualityCheck = []
    timeEnd = []

    iterationNum = 0

    print(f"Starting best x pipeline with gene name {geneName}")

    while True:

        print(f" --- Starting iteration: {iterationNum} ---")

        timeStart.append(time.time())

        shouldIncreaseX = False

        # get the new reads
        newreads, referenceReadForOrientation = getReadsIteration(
            geneName, iterationNum, referenceReadForOrientation)

        # update sample
        sampleReads += newreads
        sampleReads.sort(key=qualityReadForSorting)
        sampleReads = sampleReads[-x:]

        timeAfterGettingReads.append(time.time())

        # create consensus based on those
        consensus, worked = createNewConsensus(sampleReads)
        if not worked:
            x = int(1.5 * x)
            print(f"Increasing x to: {x}")

            iterationNum += 1

            continue

        timeAfterConsensus.append(time.time())

        # save data from this iteration
        coverages = saveQuality(
            iterationNum, allQualityConsensus, allAvgQualityReads, sampleReads)

        timeAfterQualityCheck.append(time.time())

        # get identification
        name, cov, iden = getIdentification(consensus, geneName)

        timeEnd.append(time.time())
        visualiseExecutionTime(outputDir, timeStart, timeAfterGettingReads, timeAfterConsensus, timeAfterQualityCheck, timeEnd)

        # print results to file
        outputFile.write(f"iteration: {iterationNum}\n")
        outputFile.write(f"name entry: {name}\n")
        outputFile.write(f"coverage: {cov}\n")
        outputFile.write(f"identity: {iden}\n")
        outputFile.write(f"consensus: \n{consensus}\n")
        outputFile.write("\n")
        outputFile.flush()

        # check if you can stop the pipeline
        canStop = checkEarlyStoppingCriteria(
            iterationNum, allQualityConsensus, coverages)
        if canStop:
            print("Early Stopping criterion met. Stopping pipeline")
            break

        iterationNum += 1



def main():

    geneName = sys.argv[1]
    start(geneName)

if __name__ == '__main__':
    main()