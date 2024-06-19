"""
This script is the best-x pipeline. See the README of the repo to see how to use it.
"""

from dataRetriever import *
from dataCleaner import *
from consensusMedaka import *
from identification import identify
from qualityConsensus import *
from visualise import *
import os

# minimum depth coverage that the median depth coverage along the consensus before stopping
minCoverageDepth = 40

# initial value of x
x = int(1.5 * minCoverageDepth)

# at what percent in increase in consensus quality do we stop the pipeline ?
thresholdEarlyStopping = 5

# what gene are we creating a consensus on ?
# the value must be `matK` or `rbcL` or `psbA-trnH` or `ITS`.
geneName = "matK"

# path where the new files are stored
# if it is a simulation using simulateRealTimeOutput, use "./fastqpass/" as value
readDir = "./fastqpass/"

# this variable says whether the files of the reads are generated from the simulateRealTimeOutput script
# if it is a simulation, the reads have as name the ones set by the script, and we only take 1 file per iteration
# if it isn't a simulation, we take all the new reads we can find in the folder
isSimulation = True

outputDir = "./outputPipelineBest/"


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


def getReadsIteration(geneName, iterationNum, knownFileNames, referenceReadForOrientation):

    # retrieve new reads
    if isSimulation:
        reads = waitForReadsFromFile(
            f"{readDir}iteration{iterationNum}.fastq")
        newKnownFileNames = knownFileNames
    else:
        reads, newKnownFileNames = waitForNewReadsFile(readDir, knownFileNames)

    cleanReads, referenceReadForOrientationNew = getCleanReads(
        reads, geneName, referenceReadForOrientation)

    return cleanReads, referenceReadForOrientationNew, newKnownFileNames


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

    saveCoverageVisualisation(
        outputDir, qualityConsensus, coverages, x, iterationNum)

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

    try:
        name, cov, iden = identify(pathCons, db)
        return True, name, cov, iden
    except:
        return False, "", 0, 0


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


def start():
    global x

    referenceReadForOrientation = None
    sampleReads = []

    # this is were most of the results are posted every iteration
    outputFile = open(f"{outputDir}result.txt", "w+")

    knownFileNames = set()

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

        # get the new reads
        newreads, referenceReadForOrientation, knownFileNames = getReadsIteration(
            geneName, iterationNum, knownFileNames, referenceReadForOrientation)

        # update sample
        sampleReads += newreads
        sampleReads.sort(key=qualityReadForSorting)
        sampleReads = sampleReads[-x:]

        timeAfterGettingReads.append(time.time())

        # create consensus based on those
        consensus, consensusworked = createNewConsensus(sampleReads)
        if not consensusworked:
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
        identWorked, name, cov, iden = getIdentification(consensus, geneName)

        timeEnd.append(time.time())
        visualiseExecutionTime(outputDir, timeStart, timeAfterGettingReads,
                               timeAfterConsensus, timeAfterQualityCheck, timeEnd)

        # print results to file
        outputFile.write(f"iteration: {iterationNum}\n")

        if identWorked:
            outputFile.write(f"name entry: {name}\n")
            outputFile.write(f"coverage: {cov}\n")
            outputFile.write(f"identity: {iden}\n")
        else:
            outputFile.write(f"No identification possible\n")

        outputFile.write(f"consensus: \n{consensus}\n")
        outputFile.write("\n")
        outputFile.flush()

        # check if you can stop the pipeline
        if checkEarlyStoppingCriteria(
                iterationNum, allQualityConsensus, coverages):
            print("Early Stopping criterion met! Stopping pipeline")
            break

        iterationNum += 1


def main():

    if not os.path.exists(outputDir):
        os.makedirs(outputDir)
    else:
        os.system(f"rm {outputDir}*")

    start()


if __name__ == '__main__':
    main()
