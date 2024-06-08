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
    os.system("rm ./outputMosdepth/.*")

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

    # pathAlignment = sys.argv[1]
    # quality, coverage = getQualityConsensus(pathAlignment)

    first = [13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 11, 13, 16, 16, 16, 16, 17, 17, 17, 17, 17, 17, 17, 17, 17, 17, 16, 17, 17, 17, 17, 17, 16, 17, 18, 19, 17, 15, 20, 20, 20, 20, 24, 24, 24, 24, 23, 24, 24, 24, 24, 23, 24, 22, 22, 24, 23, 24, 24, 24, 24, 24, 24, 24, 24, 26, 26, 26, 26, 25, 25, 25, 24, 24, 25, 25, 25, 25, 25, 26, 26, 26, 26, 26, 26, 25, 26, 26, 26, 25, 24, 26, 23, 22, 27, 25, 26, 25, 26, 27, 27, 28, 28, 28, 28, 28, 22, 22, 21, 26, 27, 28, 28, 28, 24, 26, 26, 28, 27, 28, 28, 28, 25, 26, 27, 28, 28, 28, 28, 28, 28, 28, 28, 25, 28, 28, 28, 27, 28, 28, 28, 27, 28, 28, 28, 28, 26, 28, 28, 27, 27, 25, 25, 28, 28, 28, 28, 27, 27, 25, 23, 24, 27, 27, 25, 26, 23, 25, 27, 27, 27, 27, 27, 26, 26, 4, 4, 25, 24, 25, 25, 25, 27, 28, 28, 25, 25, 28, 28, 27, 28, 27, 27, 28, 28, 28, 29, 29, 30, 30, 30, 30, 29, 30, 30, 30, 30, 29, 29, 26, 29, 29, 27, 29, 29, 29, 29, 30, 25, 4, 0, 1, 29, 29, 29, 29, 29, 29, 29, 28, 28, 28, 28, 28, 29, 28, 29, 30, 29, 28, 28, 31, 30, 28, 34, 34, 36, 35, 35, 36, 35, 36, 35, 35, 36, 35, 32, 36, 35, 35, 35, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 35, 35, 35, 34, 34, 35, 35, 31, 33, 30, 34, 35, 34, 37, 37, 38, 39, 39, 38, 39, 39, 39, 39, 39, 38, 37, 37, 37, 40, 38, 39, 39, 39, 39, 39, 39, 40, 40, 41, 41, 36, 36, 36, 35, 36, 36, 36, 41, 41, 36, 35, 39, 41, 41, 40, 41, 41, 41, 41, 40, 40, 41, 41, 40, 41, 41, 39, 41, 36, 37, 38, 39, 34, 36, 28, 31, 33, 38, 38, 38, 37, 38, 38, 38, 34, 28, 38, 39, 39, 39, 40, 40, 40, 38, 40, 40, 40, 39, 38, 38, 39, 39, 38, 38, 40, 40, 39, 40, 40, 34, 33, 38, 34, 37, 38, 37, 36, 36, 37, 39, 35, 37, 38, 38, 37, 38, 39, 39, 40, 40, 37, 38, 34, 37,
             38, 39, 38, 37, 36, 35, 38, 38, 33, 33, 38, 36, 37, 38, 37, 35, 34, 26, 37, 35, 37, 37, 37, 38, 38, 38, 37, 37, 37, 37, 34, 32, 33, 36, 37, 36, 37, 29, 36, 32, 24, 36, 36, 37, 37, 36, 37, 37, 33, 34, 35, 36, 36, 36, 36, 37, 37, 36, 37, 37, 36, 36, 35, 35, 35, 36, 37, 37, 36, 36, 36, 34, 30, 29, 36, 37, 37, 37, 37, 37, 37, 37, 35, 35, 34, 36, 36, 37, 36, 36, 37, 36, 37, 37, 39, 37, 38, 34, 34, 38, 38, 37, 39, 38, 39, 37, 37, 40, 40, 38, 39, 40, 40, 39, 40, 40, 41, 41, 41, 40, 40, 40, 39, 40, 39, 40, 34, 38, 39, 41, 41, 40, 41, 41, 41, 41, 41, 41, 41, 40, 40, 38, 39, 41, 41, 36, 41, 41, 40, 40, 40, 39, 40, 40, 41, 40, 40, 41, 41, 37, 40, 40, 39, 39, 39, 40, 40, 41, 41, 40, 41, 41, 41, 40, 40, 39, 34, 36, 40, 40, 40, 39, 38, 37, 37, 39, 39, 39, 39, 31, 38, 39, 37, 37, 37, 37, 35, 37, 37, 35, 35, 35, 36, 36, 36, 31, 33, 35, 36, 36, 35, 34, 34, 33, 32, 35, 34, 35, 35, 32, 32, 34, 34, 30, 25, 25, 33, 34, 34, 34, 31, 34, 34, 34, 33, 34, 34, 34, 33, 33, 33, 33, 33, 32, 33, 33, 33, 33, 33, 31, 29, 31, 32, 32, 30, 31, 33, 32, 32, 30, 32, 30, 31, 31, 31, 31, 30, 30, 30, 28, 30, 28, 29, 30, 30, 28, 30, 30, 30, 30, 29, 30, 30, 30, 29, 30, 30, 30, 19, 29, 30, 30, 30, 30, 29, 24, 28, 29, 29, 29, 29, 30, 23, 28, 28, 28, 29, 29, 29, 11, 10, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

    second = [11, 11, 11, 12, 12, 12, 12, 13, 13, 13, 13, 13, 13, 13, 13, 12, 13, 13, 13, 13, 13, 13, 13, 14, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 16, 16, 16, 16, 16, 16, 16, 16, 15, 16, 16, 16, 16, 16, 16, 15, 13, 15, 15, 14, 13, 16, 16, 13, 13, 3, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 21, 21, 21, 19, 33, 33, 22, 33, 22, 23, 26, 38, 37, 40, 38, 38, 40, 40, 40, 40, 40, 41, 41, 41, 40, 41, 41, 41, 40, 40, 40, 40, 39, 40, 40, 42, 42, 43, 43, 43, 43, 43, 43, 44, 44, 44, 44, 43, 42, 43, 42, 41, 42, 42, 43, 43, 42, 44, 44, 44, 44, 44, 44, 43, 44, 44, 44, 43, 44, 45, 42, 41, 44, 44, 44, 43, 45, 45, 45, 46, 46, 46, 45, 44, 41, 42, 41, 42, 46, 47, 48, 48, 47, 46, 46, 49, 47, 49, 49, 48, 35, 44, 47, 48, 48, 47, 48, 48, 48, 49, 49, 44, 48, 47, 48, 49, 49, 49, 49, 45, 49, 49, 48, 48, 47, 48, 48, 48, 48, 49, 48, 49, 48, 48, 47, 46, 48, 46, 46, 47, 48, 49, 48, 48, 44, 49, 50, 49, 49, 49, 50, 49, 48, 49, 43, 46, 47, 47, 46, 50, 51, 46, 46, 50, 50, 48, 49, 48, 47, 51, 51, 49, 49, 47, 52, 52, 53, 52, 52, 52, 53, 53, 53, 51, 50, 51, 47, 51, 53, 49, 52, 53, 50, 53, 52, 53, 53, 53, 54, 54, 54, 54, 52, 51, 51, 48, 48, 50, 51, 48, 47, 53, 53, 52, 52, 54, 53, 51, 54, 54, 56, 55, 54, 54, 52, 54, 54, 50, 55, 55, 44, 55, 54, 54, 55, 56, 56, 56, 54, 55, 55, 56, 55, 55, 56, 56, 56, 56, 55, 54, 56, 56, 53, 54, 49, 56, 56, 53, 57, 57, 57, 58, 58, 57, 58, 58, 57, 57, 57, 55, 57, 54, 55, 58, 56, 56, 57, 56, 57, 57, 57, 57, 57, 59, 59, 47, 47, 47, 46, 47, 47, 47, 59, 59, 46, 45, 58, 57, 57, 56, 59, 59, 59, 59, 59, 59, 59, 59, 59, 59, 59, 54, 57, 57, 58, 58, 59,
              47, 58, 43, 43, 46, 59, 59, 58, 57, 58, 58, 58, 55, 39, 58, 58, 58, 57, 58, 58, 58, 57, 59, 59, 59, 59, 58, 57, 59, 59, 58, 55, 59, 59, 58, 58, 58, 47, 46, 56, 51, 53, 58, 55, 56, 56, 56, 57, 54, 56, 58, 58, 57, 58, 58, 58, 58, 58, 57, 57, 46, 55, 56, 57, 57, 56, 55, 52, 57, 58, 53, 53, 58, 58, 58, 58, 57, 57, 55, 49, 58, 53, 56, 56, 55, 57, 57, 57, 56, 57, 57, 56, 58, 56, 58, 58, 59, 59, 59, 52, 59, 51, 44, 59, 57, 59, 57, 57, 59, 59, 53, 53, 56, 58, 58, 59, 59, 58, 59, 59, 58, 59, 58, 58, 56, 57, 58, 57, 57, 57, 57, 58, 57, 53, 40, 41, 57, 58, 57, 57, 57, 56, 56, 56, 56, 54, 53, 54, 54, 53, 53, 53, 53, 52, 53, 53, 52, 47, 50, 48, 47, 51, 52, 51, 51, 51, 51, 50, 50, 53, 53, 48, 50, 51, 51, 52, 50, 51, 52, 53, 53, 51, 50, 52, 52, 52, 51, 52, 48, 52, 48, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 51, 50, 50, 51, 52, 51, 49, 51, 52, 51, 51, 51, 49, 51, 51, 52, 52, 51, 52, 52, 47, 52, 52, 49, 51, 51, 51, 51, 51, 51, 49, 50, 49, 49, 50, 51, 50, 47, 49, 51, 51, 51, 50, 49, 49, 49, 50, 50, 49, 49, 35, 50, 49, 46, 46, 45, 45, 45, 44, 44, 45, 45, 44, 44, 44, 44, 38, 43, 44, 43, 44, 44, 41, 43, 42, 41, 44, 44, 44, 44, 43, 43, 43, 43, 39, 35, 32, 43, 43, 42, 41, 38, 42, 42, 42, 41, 41, 39, 39, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 38, 35, 35, 38, 37, 36, 38, 38, 38, 37, 36, 37, 36, 37, 37, 37, 37, 37, 36, 36, 34, 36, 36, 36, 36, 35, 36, 36, 36, 36, 35, 36, 36, 36, 36, 36, 36, 36, 35, 16, 33, 34, 34, 34, 34, 34, 29, 33, 33, 32, 33, 32, 33, 21, 33, 32, 33, 33, 31, 31, 28, 16, 13, 8, 3, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

    third = [28, 28, 28, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 28, 28, 29, 28, 28, 25, 27, 28, 28, 26, 28, 29, 31, 35, 35, 35, 35, 36, 36, 36, 36, 36, 35, 37, 37, 38, 38, 39, 39, 39, 39, 39, 39, 39, 39, 38, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 37, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 38, 38, 38, 39, 39, 39, 39, 36, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 38, 39, 38, 38, 38, 38, 38, 35, 36, 51, 26, 27, 27, 50, 50, 50, 50, 50, 50, 51, 51, 50, 50, 52, 53, 51, 52, 51, 52, 53, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 53, 53, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 51, 54, 54, 54, 54, 54, 55, 54, 53, 55, 54, 54, 55, 55, 55, 55, 56, 56, 56, 56, 54, 54, 55, 55, 49, 28, 55, 57, 58, 58, 58, 57, 57, 58, 57, 58, 58, 58, 48, 55, 56, 56, 56, 56, 57, 57, 57, 58, 58, 54, 57, 56, 57, 58, 58, 58, 58, 54, 58, 58, 57, 57, 53, 57, 57, 56, 56, 57, 56, 57, 39, 39, 39, 49, 52, 53, 50, 52, 53, 56, 56, 55, 57, 58, 57, 54, 56, 57, 56, 57, 58, 58, 58, 57, 58, 54, 56, 56, 56, 55, 56, 58, 54, 56, 58, 58, 53, 58, 57, 57, 58, 58, 58, 56, 55, 57, 57, 58, 58, 58, 57, 58, 58, 58, 57, 58, 59, 56, 57, 59, 56, 58, 59, 58, 58, 57, 58, 58, 58, 58, 58, 59, 59, 57, 58, 58, 56, 57, 57, 57, 56, 56, 59, 59, 58, 58, 59, 59, 59, 58, 58, 59, 59, 59, 59, 59, 59, 59, 58, 59, 59, 35, 59, 58, 57, 59, 59, 59, 59, 57, 57, 59, 59, 59, 59, 59, 59, 59, 59, 58, 58, 59, 59, 57, 57, 55, 59, 59, 57, 59, 57, 57, 58, 59, 58, 59, 59, 59, 59, 59, 57, 58, 55, 58, 59, 59, 59, 59, 59, 59, 58, 58, 59, 59, 59, 59, 35, 35, 35, 34, 35, 35, 35, 59, 59, 35, 33, 57,
             58, 59, 56, 59, 59, 59, 59, 59, 59, 59, 59, 59, 59, 59, 56, 59, 57, 58, 58, 59, 34, 58, 33, 34, 35, 59, 58, 58, 57, 59, 58, 59, 59, 31, 59, 58, 58, 58, 58, 58, 58, 58, 59, 59, 59, 59, 59, 58, 58, 59, 59, 58, 59, 59, 58, 59, 59, 35, 34, 56, 53, 56, 58, 57, 59, 58, 58, 59, 54, 57, 59, 59, 58, 59, 59, 59, 59, 59, 59, 59, 35, 59, 58, 58, 58, 58, 57, 58, 59, 59, 56, 55, 58, 59, 59, 58, 58, 57, 58, 54, 57, 53, 56, 56, 55, 58, 58, 58, 58, 57, 56, 58, 58, 55, 59, 59, 59, 59, 58, 55, 59, 57, 31, 59, 56, 58, 56, 57, 59, 59, 55, 55, 59, 59, 58, 59, 59, 59, 59, 59, 59, 59, 59, 59, 58, 58, 59, 59, 59, 59, 59, 59, 59, 57, 59, 59, 59, 59, 59, 59, 59, 59, 59, 59, 59, 59, 58, 58, 58, 58, 58, 58, 58, 57, 57, 55, 56, 54, 54, 58, 58, 58, 58, 57, 58, 57, 57, 57, 58, 54, 54, 56, 57, 58, 57, 57, 58, 58, 58, 58, 57, 57, 57, 57, 57, 58, 57, 56, 54, 56, 57, 58, 57, 57, 58, 58, 58, 57, 57, 58, 57, 56, 56, 57, 55, 55, 57, 57, 57, 57, 57, 55, 56, 54, 55, 57, 57, 57, 56, 31, 53, 56, 53, 53, 54, 54, 55, 55, 55, 55, 54, 55, 54, 54, 53, 53, 54, 51, 51, 55, 55, 55, 55, 55, 54, 54, 54, 54, 53, 52, 26, 52, 52, 48, 48, 48, 48, 48, 48, 48, 47, 48, 47, 47, 48, 48, 42, 46, 48, 47, 47, 48, 46, 48, 47, 47, 48, 48, 48, 48, 47, 47, 46, 47, 44, 41, 39, 47, 47, 47, 46, 44, 46, 46, 46, 46, 46, 45, 45, 44, 44, 44, 43, 43, 42, 43, 42, 42, 42, 42, 42, 38, 38, 42, 41, 39, 42, 42, 42, 41, 42, 42, 42, 42, 42, 42, 41, 41, 40, 40, 39, 39, 39, 39, 39, 39, 39, 38, 38, 38, 38, 38, 26, 26, 24, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 36, 33, 36, 36, 36, 36, 36, 36, 36, 34, 34, 34, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

    saveCoverageVisualisationTOGETHER(first, second, third)


if __name__ == "__main__":
    main()
