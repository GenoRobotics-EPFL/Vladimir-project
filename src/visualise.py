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


def visualiseBaseQualities(pathFastqFile, readQualityConverted):

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

    visualiseBaseQualities(pathFastqFile, readQualityConverted)

def visualiseReadQualities(pathFastqFile, readQualities):

    plt.figure(figsize=(10, 6))

    plt.hist(readQualities, bins=range(0, 95))

    # plot the mean
    mean = sum(readQualities) / len(readQualities)
    plt.axvline(mean, color="r", label="mean :" + str(round(mean, 2)))

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

    readQualities = [ getQualityRead(read) for read in reads]

    visualiseReadQualities(pathFastqFile, readQualities)





def main():

    pathFastqFile = sys.argv[1]

    # baseQualityDistribution(pathFastqFile)

    readQualityDistribution(pathFastqFile)


if __name__ == "__main__":
    main()
