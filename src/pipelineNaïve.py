import sys
from dataRetriever import *
from dataCleaner import *
from consensusMedaka import *
from identification import identify
from qualityConsensus import *
from visualise import *

# used for preprocessing
referenceReadForOrientation = None
allReads = []

# at what percent in increase in consensus quality do we stop the pipeline ?
thresholdEarlyStopping = 5

outputDir = "./outputPipelineNaive/"
os.system(f"rm {outputDir}*")  # clean dir to make sure we start fresh

# this is were most of the results are posted every iteration
outputFile = open(f"{outputDir}results.txt", "w+")

# this is used to plot the quality of the reads in the sample over time
allAvgQualityReads = []
# this is used to plot the quality of the consensus over time
allQualityConsensus = []

timeStart = []
timeAfterGettingReads = []
timeAfterConsensus = []
timeAfterQualityCheck = []
timeEnd = []


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


def getReadsLaterIterations(geneName, iterationNum):

    # retrieve new reads
    reads = waitForReadsFromFile(f"fastqpass/iteration{iterationNum}.fastq")

    cleanReads, _ = getCleanReads(
        reads, geneName, referenceReadForOrientation)

    # append new reads at the front as they are more likely to have a low quality
    global allReads
    allReads += cleanReads


def splitReads2080():

    allReads.sort(key=qualityReadForSorting)

    splitIndex = int(len(allReads) * 0.2)

    best20Reads = allReads[-splitIndex:]
    worse80Reads = allReads[:-splitIndex]

    return best20Reads, worse80Reads


def createNewConsensus(best20, worse80):

    # write those reads to file so consensus can act on those reads
    pathBest20 = f"{outputDir}temp_reads20ForConsensus.fasta"
    writeReadsToFile(pathBest20, best20)

    pathallReads = f"{outputDir}temp_readsAllForConsensus.fasta"
    writeReadsToFile(pathallReads, allReads)

    consensus = getConsensus8020(pathBest20, pathallReads)

    return consensus


def saveQuality(iterationNum):

    qualityConsensus, coverages = getQualityConsensus(
        "./outputMedaka/calls_to_draft.bam")
    allQualityConsensus.append(qualityConsensus)

    saveCoverageVisualisation(outputDir,
                              qualityConsensus, coverages, len(allReads), iterationNum)

    qualityReads = sum(map(qualityReadForSorting, allReads)) / len(allReads)
    allAvgQualityReads.append(qualityReads)

    saveQualityReadsOverTime(outputDir, allAvgQualityReads)
    saveQualityConsensusOverTime(outputDir, allQualityConsensus)
    saveQualityBothOverTime(outputDir, allAvgQualityReads, allQualityConsensus)

    return coverages


def getIdentification(consensus, db):

    pathCons = f"{outputDir}tempConsForIdent.fasta"
    writeConsensus(pathCons, consensus)

    # db = "rbcL" # "psbA-trnH"
    name, cov, iden = identify(pathCons, db)

    return name, cov, iden


def checkEarlyStoppingCriteria(iterationNum, coverages):

    secondCondition = np.median(coverages) >= 40

    return secondCondition


def start(geneName):

    print(f"Starting naive pipeline with gene name {geneName}")

    iterationNum = 0

    while True:

        print(f" --- Starting iteration: {iterationNum} ---")
        timeStart.append(time.time())

        # get the new reads
        if iterationNum == 0:
            getReadsFirstIteration(geneName)
        else:
            getReadsLaterIterations(geneName, iterationNum)

        timeAfterGettingReads.append(time.time())

        best20, worse80 = splitReads2080()

        # create consensus based on those
        consensus = createNewConsensus(best20, worse80)

        timeAfterConsensus.append(time.time())

        # save data from this iteratoin
        coverages = saveQuality(iterationNum)

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
        canStop = checkEarlyStoppingCriteria(iterationNum, coverages)
        if canStop:
            print("Early Stopping criterion met. Stopping pipeline")
            break

        iterationNum += 1


def main():

    geneName = sys.argv[1]
    start(geneName)


if __name__ == "__main__":
    main()
