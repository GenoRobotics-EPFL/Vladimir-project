"""
This script is used to calculate the quality of a consensus.
"""

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

    plt.figure()
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
    plt.axhline(y=40, color="r", label="threshold coverage: 40")

    # plot the quality information of consensus using metric
    start = 0
    for i, cov in enumerate(coverages):
        if cov >= median:
            start = i
            break

    lenConsAbove = lengthConsensusAboveMedian(coverages, median)
    plt.plot([start, start + lenConsAbove],
             [median, median], c="g", label=f"portion above median with length: {lenConsAbove}")
    # plt.plot([], [], ' ', label=f"quality consensus: {qualityConsensus}")

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
    pathToOutputDir = "./outputMosdepth"
    if not os.path.exists(pathToOutputDir):
        os.makedirs(pathToOutputDir)
    else:
        os.system(f"rm {pathToOutputDir}/*")

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


def saveCoverageVisualisationTOGETHER(coverage1, coverage2, coverage3):

    plt.figure()
    plt.title(f"Evolution of depth coverage of consensus:")

    coverage1 = [cov for cov in coverage1 if cov != 0]
    coverage2 = [cov for cov in coverage2 if cov != 0]
    coverage3 = [cov for cov in coverage3 if cov != 0]

    # center all of them
    maxLen = max(max(len(coverage1), len(coverage2)), len(coverage3))

    shift = int((maxLen - len(coverage3)) / 2)
    xCoverage3 = [i + shift for i in range(len(coverage3))]
    plt.plot(xCoverage3, coverage3, c="g", alpha=0.5,
             label="coverage iteration 24")

    shift = int((maxLen - len(coverage2)) / 2)
    xCoverage2 = [i + shift for i in range(len(coverage2))]
    plt.plot(xCoverage2, coverage2, c="tab:orange", alpha=0.5,
             label="coverage iteration 10")

    shift = int((maxLen - len(coverage1)) / 2)
    xCoverage1 = [i + shift for i in range(len(coverage1))]
    plt.plot(xCoverage1, coverage1, c="b",
             alpha=0.5, label="coverage iteration 3")

    # # plot the median
    # median1 = np.median(coverage1)
    # plt.axhline(y=median1, color="b")

    # median2 = np.median(coverage2)
    # plt.axhline(y=median2, alpha=0.2, color="tab:orange")

    # median3 = np.median(coverage3)
    # plt.axhline(y=median3, color="g")

    # plot the min we allow
    plt.axhline(y=40, color="r", label="minDepthCoverage: 40")

    plt.legend(loc="lower left")

    plt.ylabel("depth coverage")
    plt.xlabel("base in consensus")

    plt.savefig(f"coverageIterationTOGETHER.png")


def main():

    pathAlignment = sys.argv[1]
    quality, coverage = getQualityConsensus(pathAlignment)

    print(quality)

if __name__ == "__main__":
    main()
