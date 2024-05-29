import sys
import os
from dataRetriever import *
import matplotlib.pyplot as plt
import numpy as np


def getCoverageFromFile(pathBedFile):

    lines = getLinesFromFile(pathBedFile)

    coverageReads = []

    for line in lines:

        values = line.split("\t")

        start = int(values[1])
        end = int(values[2])
        coverage = int(values[3])

        for _ in range(end - start):
            coverageReads.append(coverage)

    return coverageReads


def visualiseCoverage(baseCoverages):

    plt.figure(figsize=(10, 6))

    plt.plot(baseCoverages)

    plt.ylabel("Coverage")

    plt.xlabel("base on consensus")

    plt.show()


def saveCoverageVisualisation(outputDir, qualityConsensus, coverages, numReadsInConsensus, iterationNum):

    plt.figure(figsize=(10, 6))
    plt.title(f"Coverage of consensus for iteration: {iterationNum}")

    plt.plot(coverages)

    # show how many reads were used
    plt.ylim(0, numReadsInConsensus + 1)
    plt.axhline(y=numReadsInConsensus, color="k",
                label=f"x :{numReadsInConsensus}")

    # plot the median
    median = np.median(coverages)
    plt.axhline(y=median, color="b",
                label=f"median : {round(median, 2)}")

    # plot the min we allow
    plt.axhline(y=40, color="r")

    # plot the quality information of consensus using metric
    start = 0
    for i, cov in enumerate(coverages):
        if cov >= median:
            start = i
            break

    lenConsAbove = lengthConsensusAboveMedian(coverages, median)
    plt.plot([start, start + lenConsAbove],
             [median, median], c="g", label=f"portion above median with length: {lenConsAbove}")
    plt.plot([], [], ' ', label=f"quality consensus: {qualityConsensus}")

    plt.legend(loc="lower left")

    plt.ylabel("depth coverage")
    plt.xlabel("base in consensus")

    plt.savefig(f"{outputDir}coverageIteration{iterationNum}.png")


def getValidityCoverage(coverages, numReads):

    coverages.sort()

    portionLowest = coverages[:int(len(coverages) * .3)]
    avgLowest = sum(portionLowest) / len(portionLowest)

    threshold = 0.3 * numReads

    return avgLowest > threshold and avgLowest > 5


def getCoverageFromConsensus(pathAlignment):

    # first run mosdepth
    command = f"./src/mosdepth ./outputMosdepth/ {pathAlignment}"
    os.system(command)

    # uncompress output data from mosdepth
    command = "gzip -f -dk ./outputMosdepth/.per-base.bed.gz"
    os.system(command)

    coverages = getCoverageFromFile("./outputMosdepth/.per-base.bed")

    return coverages


def lengthConsensusAboveMedian(coverages, medianCoverage):

    # find first base with coverage above median
    start = 0
    for i, cov in enumerate(coverages):
        if cov >= medianCoverage:
            start = i
            break

    # find last base with coverage above median
    end = 0
    for i, cov in reversed(list(enumerate(coverages))):
        if cov >= medianCoverage:
            end = i
            break

    return end - start


def getQualityConsensus(pathAlignment):

    coverages = getCoverageFromConsensus(pathAlignment)

    medianCoverage = np.median(coverages)

    lenConsensus = lengthConsensusAboveMedian(coverages, medianCoverage)

    return medianCoverage * lenConsensus, coverages


def main():

    pathAlignment = sys.argv[1]

    quality, coverage = getQualityConsensus(pathAlignment)

    print(quality)


if __name__ == "__main__":
    main()
