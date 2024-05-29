import sys
import os
import matplotlib.pyplot as plt
from dataRetriever import *
from dataCleaner import *


def getAllQualities(reads):

    qualities = []

    for read in reads:
        for qualityCh in read.quality:

            qualities.append(qualityCh)

    return qualities


def convertQualities(allReadsQuality):

    # note that here we can't distinguish the reads from each other anymore
    # they are all in one list
    convertedQualities = []

    for qualityCh in allReadsQuality:

        # transform the character into the int quality
        score = ord(qualityCh) - 33

        convertedQualities.append(score)

    return convertedQualities


def visualiseBaseQualityDitribution(pathFastqFile, readQualityConverted):

    plt.figure(figsize=(10, 6))

    plt.hist(readQualityConverted, bins=range(0, 95))

    # plot the mean
    mean = sum(readQualityConverted) / len(readQualityConverted)
    plt.axvline(mean, color="r", label="mean :" + str(round(mean, 2)))

    plt.legend()

    plt.ylabel("Frequency")

    plt.xlabel("Base quality score 0-93 (higher is better)")
    plt.xlim(0, 95)

    nameFile = os.path.basename(pathFastqFile)
    plt.title(nameFile)

    nameOutput = nameFile + ".png"
    plt.savefig(nameOutput)  # saves in current directory


def baseQualityDistribution(pathFastqFile):

    reads = getReadsFromFile(pathFastqFile)

    allQualities = getAllQualities(reads)

    readQualityConverted = convertQualities(allQualities)

    visualiseBaseQualityDitribution(pathFastqFile, readQualityConverted)


def visualiseReadQualityDistribution(pathFastqFile, readQualities):

    plt.figure(figsize=(10, 6))

    plt.hist(readQualities, bins=range(0, 95))

    # plot the mean
    mean = sum(readQualities) / len(readQualities)
    plt.axvline(mean, color="b", label="mean :" + str(round(mean, 2)))

    # plot the mim value tha tis removed when preprocessing
    plt.axvline(10, color="r", label="min : 10")

    plt.legend()

    plt.ylabel("Frequency")

    plt.xlabel("Read quality score 0-93 (higher is better)")
    plt.xlim(0, 95)

    nameFile = os.path.basename(pathFastqFile)
    plt.title(nameFile)

    nameOutput = nameFile + ".png"
    plt.savefig(nameOutput)  # saves in current directory


def readQualityDistribution(pathFastqFile):

    reads = getReadsFromFile(pathFastqFile)

    readQualities = [getQualityRead(read) for read in reads]

    visualiseReadQualityDistribution(pathFastqFile, readQualities)


def visualiseReadLengthDistribution(pathFastqFile, readLengths):

    plt.figure(figsize=(10, 6))

    plt.hist(readLengths, bins=30)

    # mean
    mean = sum(readLengths) / len(readLengths)
    plt.axvline(mean, color="b", label=f"mean : {mean}")

    plt.axvline(300, color="r", label=f"min : 300")
    plt.axvline(850, color="r", label=f"max : 850")

    plt.ylabel("Frequency")
    plt.xlabel("Read length")

    plt.legend()

    nameFile = os.path.basename(pathFastqFile)
    plt.title(nameFile)

    plt.savefig(nameFile + "lengthDist.png")


def readLengthDistribution(pathFastqFile):

    reads = getReadsFromFile(pathFastqFile)

    readLengths = [len(read.sequence) for read in reads]

    visualiseReadLengthDistribution(pathFastqFile, readLengths)


def saveQualityConsensusOverTime(outputDir, allqualityconsensus):

    plt.figure(figsize=(10, 6))
    plt.title(f"Evolution of consensus quality over time:")
    plt.plot(allqualityconsensus)
    plt.xlabel("Iteration")
    plt.ylabel("consensus quality")

    # plot the early sotpping line
    if (len(allqualityconsensus) >= 5):

        x1 = len(allqualityconsensus) - 5
        x2 = len(allqualityconsensus) - 1
        y1 = allqualityconsensus[x1]
        y2 = allqualityconsensus[x2]

        maxQual = max(allqualityconsensus)
        minQual = min(allqualityconsensus)

        percentIncreaseSinceStart = (y2 - y1) / (maxQual - minQual)

        plt.plot([x1, x2], [y1, y2], c="r",
                 label=f"early stopping check: {int(percentIncreaseSinceStart * 100)}%")

        plt.legend()

    plt.savefig(f"{outputDir}evolutionConsensusQuality.png")
    plt.close()


def saveQualityReadsOverTime(outputDir, allqualityreads):

    # first save quality of read over time
    plt.figure(figsize=(10, 6))
    plt.title(f"Evolution of average read quality over time")
    plt.plot(allqualityreads)
    plt.xlabel("Iteration")
    plt.ylabel("read quality")
    plt.savefig(f"{outputDir}evolutionReadQuality.png")
    plt.close()


def visualiseExecutionTime(outputDir, timeStart, timeAfterGettingReads,
                           timeAfterConsensus, timeAfterQualityCheck, timeEnd):

    # total time per iteration
    timeTotal = []
    for i in range(len(timeStart)):
        timeTotal.append(timeEnd[i] - timeStart[i])

    plt.figure(figsize=(10, 6))
    plt.title(f"Time taken for each step:")
    plt.bar(list(range(len(timeStart))), timeTotal,
            color="y", label="time identification")
    plt.xlabel("Iteration")
    plt.ylabel("time (s)")
    plt.savefig(f"{outputDir}executionTime.png")
    plt.close()


def saveQualityBothOverTime(outputDir, allqualityreads, allqualityconsensus):

    # now draw quality of consensus over time as quality of read improves
    fig, ax1 = plt.subplots()
    color = 'tab:red'
    ax1.set_xlabel('iteration number')
    ax1.set_ylabel('average quality of reads in sample', color=color)
    ax1.plot(allqualityreads, color=color)
    ax1.tick_params(axis='y', labelcolor=color)

    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

    # we already handled the x-label with ax1
    color = 'tab:blue'
    ax2.set_ylabel('quality of consensus', color=color)
    ax2.plot(allqualityconsensus, color=color)
    ax2.tick_params(axis='y', labelcolor=color)
    fig.tight_layout()  # otherwise the right y-label is slightly clipped

    plt.savefig(f"{outputDir}evolutionQuality.png")
    plt.close()


def main():

    pathFastqFile = sys.argv[1]

    readLengthDistribution(pathFastqFile)
    readQualityDistribution(pathFastqFile)


if __name__ == "__main__":
    main()
