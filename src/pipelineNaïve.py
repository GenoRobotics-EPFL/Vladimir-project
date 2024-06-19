"""
This script is the naive pipeline. See the README of the repo to see how to use it.
"""

import sys
from dataRetriever import *
from dataCleaner import *
from consensusMedaka import *
from identification import identify
from qualityConsensus import *
from visualise import *


# ! START variables to change depending on use case

# at what percent in increase in consensus quality do we stop the pipeline ?
thresholdEarlyStopping = 5

# what gene are we creating a consensus on ?
# the value must be `matK` or `rbcL` or `psbA-trnH` or `ITS`.
geneName = "ITS"

# this variable says whether the files of the reads are generated from the simulateRealTimeOutput script
# if it is a simulation, the reads have as name the ones set by the script, and we only take 1 file per iteration
# if it isn't a simulation, we take all the new reads we can find in the folder
isSimulation = True

# path where the new files are stored
# if it is a simulation using simulateRealTimeOutput, use "./fastqpass/" as value
readDir = "./fastqpass/"

# this is where all the main results are outputed
outputDir = "./outputPipelineNaive/"

# ! END variables to change depending on use case


def qualityReadFor2080(read):
    """
    This method defines what reads will be considered 'best' when performing 20 80
    """
    return len(read.sequence)


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


def splitReads2080(allReads):

    allReads.sort(key=qualityReadFor2080)

    splitIndex = int(len(allReads) * 0.2)

    best20Reads = allReads[-splitIndex:]
    worse80Reads = allReads[:-splitIndex]

    return best20Reads, worse80Reads


def createNewConsensus(best20, worse80):

    # write those reads to file so consensus can act on those reads
    pathBest20 = f"./temp_reads20ForConsensus.fasta"
    writeReadsToFile(pathBest20, best20)

    pathallReads = f"./temp_readsAllForConsensus.fasta"
    writeReadsToFile(pathallReads, best20 + worse80)

    try:
        consensus = getConsensus8020(pathBest20, pathallReads)
    except NoConsensusEx:
        return "", False

    return consensus, True


def saveQuality(iterationNum, allQualityConsensus, allAvgQualityReads, allReads):

    qualityConsensus, coverages = getQualityConsensus(
        "./outputMedaka/calls_to_draft.bam")
    allQualityConsensus.append(qualityConsensus)

    saveCoverageVisualisation(outputDir,
                              qualityConsensus, coverages, len(allReads), iterationNum)

    qualityReads = sum(map(qualityReadFor2080, allReads)) / len(allReads)
    allAvgQualityReads.append(qualityReads)

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


def checkEarlyStoppingCriteria(iterationNum, coverages):

    secondCondition = np.median(coverages) >= 40

    return secondCondition


def start():

    referenceReadForOrientation = None

    allReads = []

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

    print(f"Starting naive pipeline with gene name {geneName}")

    while True:

        print(f" --- Starting iteration: {iterationNum} ---")
        timeStart.append(time.time())

        # get the new reads
        newreads, referenceReadForOrientation, knownFileNames = getReadsIteration(
            geneName, iterationNum, knownFileNames, referenceReadForOrientation)

        allReads += newreads

        timeAfterGettingReads.append(time.time())

        # get best 20
        best20, worse80 = splitReads2080(allReads)

        # create consensus based on those
        consensus, consworked = createNewConsensus(best20, worse80)
        if not consworked:
            print(f"couldn't create consensus !")
            iterationNum += 1
            continue

        timeAfterConsensus.append(time.time())

        # save data from this iteratoin
        coverages = saveQuality(
            iterationNum, allQualityConsensus, allAvgQualityReads, allReads)

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
        canStop = checkEarlyStoppingCriteria(iterationNum, coverages)
        if canStop:
            print("Early Stopping criterion met. Stopping pipeline")
            break

        iterationNum += 1


def main():

    if not os.path.exists(outputDir):
        os.makedirs(outputDir)
    else:
        os.system(f"rm {outputDir}*")

    start()


if __name__ == "__main__":
    main()
