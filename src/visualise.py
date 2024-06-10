"""
This file contains all the graphs that were used for the report and debugging.
"""

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

    plt.figure()

    plt.hist(readQualities, bins=range(0, 95))

    # plot the mean
    mean = sum(readQualities) / len(readQualities)
    plt.axvline(mean, color="b", label="mean :" + str(round(mean, 2)))

    # plot the mim value tha tis removed when preprocessing
    plt.axvline(10, color="r", label="min : 10")

    plt.legend()

    plt.ylabel("Frequency")

    plt.xlabel("Read quality score 0-93")
    plt.xlim(0, 95)

    nameFile = os.path.basename(pathFastqFile)
    plt.title(nameFile)

    nameOutput = nameFile + ".png"
    plt.savefig(nameOutput)  # saves in current directory


def readQualityDistribution(pathFastqFile):

    reads = getReadsFromFile(pathFastqFile)

    readQualities = [getQualityRead(read) for read in reads]

    visualiseReadQualityDistribution(pathFastqFile, readQualities)


def visualiseReadLengthDistribution(pathFastqFile, readLengths, geneName):
    if geneName == "matK":
        minReadLength = 300
        avg = 700
        maxReadLength = 900

    elif geneName == "rbcL":
        minReadLength = 300
        avg = 750
        maxReadLength = 850

    elif geneName == "psbA-trnH":
        minReadLength = 300
        avg = 500
        maxReadLength = 700

    elif geneName == "ITS":
        minReadLength = 300
        avg = 700
        maxReadLength = 800

    plt.figure()

    plt.hist(readLengths, bins=30)

    # mean
    mean = sum(readLengths) / len(readLengths)
    plt.axvline(mean, color="b", label=f"mean : {round(mean, 2)}")

    plt.axvline(minReadLength, color="r", label=f"min : {minReadLength}")
    plt.axvline(maxReadLength, color="r", label=f"max : {maxReadLength}")

    plt.ylabel("Frequency")
    plt.xlabel("Read length")

    plt.legend()

    nameFile = os.path.basename(pathFastqFile)
    plt.title(nameFile)

    plt.savefig(nameFile + "lengthDist.png")


def readLengthDistribution(pathFastqFile, geneName):

    reads = getReadsFromFile(pathFastqFile)

    readLengths = [len(read.sequence) for read in reads]

    visualiseReadLengthDistribution(pathFastqFile, readLengths, geneName)


def saveQualityConsensusOverTime(outputDir, allqualityconsensus):

    plt.figure(figsize=(10, 6))
    plt.title(f"Evolution of consensus quality over time:")
    plt.plot(allqualityconsensus)
    plt.xlabel("Iteration")
    plt.ylabel("consensus quality")

    # plot the early sotpping line
    if (len(allqualityconsensus) >= 8):

        x1 = len(allqualityconsensus) - 1 - 7
        x2 = len(allqualityconsensus) - 1 - 6
        meanbefore = (allqualityconsensus[x1] + allqualityconsensus[x2]) / 2

        x1 = len(allqualityconsensus) - 1
        x2 = len(allqualityconsensus) - 2
        meanafter = (allqualityconsensus[x1] + allqualityconsensus[x2]) / 2

        maxQual = max(allqualityconsensus)
        minQual = min(allqualityconsensus)

        percentIncreaseSinceStart = (
            meanafter - meanbefore) / (maxQual - minQual)

        plt.plot([len(allqualityconsensus) - 8.5, len(allqualityconsensus) - 1.5], [meanbefore, meanafter], c="r",
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


def visualiseAvgReadQualityTogether():

    plt.figure()
    plt.title(f"Evolution of average read quality in sample over time")

    first = [109203.0, 252564.0, 395858.0, 508983.0, 606109.0, 651563.0, 677443.0, 695595.0, 707243.0, 728571.0, 736996.0, 741847.0,
             745117.0, 749345.0, 762349.0, 767330.0, 774678.0, 778177.0, 782433.0, 788155.0, 789132.0, 791743.0, 793872.0, 796831.0, 797905.0]
    plt.plot(first, label="Allium Ursinum (good)")

    second = [80065.0, 174836.0, 260289.0, 363163.0, 450814.0, 502545.0, 524965.0, 538541.0, 552130.0, 568405.0, 579409.0, 596760.0, 602262.0, 614349.0,
              623745.0, 629414.0, 634220.0, 637952.0, 647820.0, 649530.0, 652354.0, 653977.0, 657098.0, 662692.0, 665088.0, 667378.0, 669088.0, 671020.0, 671840.0]

    plt.plot(second, label="Ficus religiosa (good)")
    third = [48034.0, 119641.0, 151012.0, 182560.0, 217861.0, 270281.0, 308638.0, 339978.0, 397344.0, 416435.0, 450554.0, 505589.0, 531273.0,
             560090.0, 578481.0, 583057.0, 594870.0, 612476.0, 625940.0, 635521.0, 640455.0, 651238.0, 657050.0, 660365.0, 665044.0, 671197.0]
    plt.plot(third, label="Solanum Lycopersicum (good)")

    first = [90300.0, 203446.0, 308543.0, 395547.0, 487050.0, 529636.0, 554676.0, 567286.0, 578372.0, 598121.0, 604606.0, 608832.0, 611754.0,
             615666.0, 628384.0, 632795.0, 640410.0, 643724.0, 647233.0, 652410.0, 653187.0, 655775.0, 657854.0, 660609.0, 661699.0, 662004.0]
    plt.plot(first, label="Allium Ursinum (bad)")

    second = [51050.0, 118947.0, 172410.0, 244425.0, 307787.0, 374053.0, 405889.0, 418990.0, 427706.0, 441193.0, 447202.0, 463138.0, 468334.0, 481522.0, 489078.0,
              494938.0, 500095.0, 506050.0, 515507.0, 517682.0, 520710.0, 522222.0, 523345.0, 528448.0, 530570.0, 531500.0, 533236.0, 535442.0, 536177.0, 537242.0, 537844.0]
    plt.plot(second, label="Ficus religiosa (bad)")
    third = [30296.0, 87375.0, 105867.0, 128570.0, 154485.0, 195734.0, 224591.0, 257249.0, 302287.0, 315097.0, 340090.0, 383226.0, 401048.0, 427798.0, 440908.0, 441804.0, 446880.0, 458783.0,
             468925.0, 477469.0, 482658.0, 492504.0, 501041.0, 504872.0, 509787.0, 520227.0, 523890.0, 527462.0, 527462.0, 535807.0, 540547.0, 541370.0, 541370.0, 543534.0, 543534.0, 549568.0]
    plt.plot(third, label="Solanum Lycopersicum (bad)")

    plt.xlabel("Iteration number")
    plt.ylabel("average read quality in sample")
    plt.legend()
    plt.savefig("evolutionReadQualityALL.png")
    plt.close()


def visualiseCOnsQualityTogether():

    plt.figure()
    plt.title(f"Evolution of consensus quality over time")

    first = [3255.0, 6015.0, 13041.0, 15240.0, 14040.0, 16560.0, 17028.0, 22066.0, 17820.0, 19680.0, 23040.0, 23950.0,
             24000.0, 24990.0, 25168.0, 25228.0, 25056.0, 25920.0, 29268.0, 29645.0, 29865.0, 29975.0, 27170.0, 27225.0, 28325.0]
    plt.plot(first, label="Allium Ursinum (good)")
    second = [3210.0, 4368.0, 6006.0, 8064.0, 9622.0, 11232.0, 12888.0, 13542.0, 14022.0, 15129.0, 16154.0, 14265.0, 14265.0, 15839.0,
              15792.0, 16121.0, 16575.0, 16219.0, 17013.0, 16907.0, 17316.0, 17316.0, 17545.0, 17528.0, 17808.0, 17841.0, 22220.0, 18240.0, 18240.0]
    plt.plot(second, label="Ficus religiosa (good)")
    third = [1636.0, 3344.0, 4004.0, 5330.0, 8384.0, 8721.0, 8694.0, 8740.0, 11492.0, 10780.0, 13392.0, 14070.0, 15696.0,
             15542.0, 15918.0, 18480.0, 17415.0, 16500.0, 17955.0, 18095.0, 18473.0, 19850.0, 19656.0, 19552.0, 17649.0, 19344.0]
    plt.plot(third, label="Solanum Lycopersicum (good)")

    first = [2604.0, 4064.0, 10269.0, 14084.0, 8109.0, 15947.0, 17343.0, 16674.0, 19110.0, 18630.0, 21574.0, 22607.0, 22560.0,
             23088.0, 22800.0, 23500.0, 23460.0, 23088.0, 24327.0, 25380.0, 26784.0, 26676.0, 24894.0, 24138.0, 26244.0, 23052.0]
    plt.plot(first, label="Allium Ursinum (bad)")
    second = [1452.0, 3608.0, 5265.0, 6573.0, 8181.0, 10200.0, 10116.0, 10812.0, 11351.5, 12402.0, 12402.0, 13160.0, 12997.0, 13815.0, 14899.0,
              15134.0, 15600.0, 15552.0, 16116.0, 15242.5, 16050.0, 16371.0, 16484.0, 16744.0, 17325.0, 17160.0, 17160.0, 17545.0, 17435.0, 17752.0, 17670.0]
    plt.plot(second, label="Ficus religiosa (bad)")
    third = [1046.0, 2863.0, 2952.0, 3717.0, 2862.0, 4494.0, 4650.0, 5072.0, 4459.0, 5310.0, 0.0, 8075.0, 6946.0, 10571.0, 10890.0, 10140.0, 9833.5, 9108.0, 12312.0,
             12168.0, 13860.0, 13860.0, 13690.0, 13690.0, 14060.0, 15048.0, 14858.0, 14937.0, 14937.0, 14880.0, 15352.0, 14781.0, 14781.0, 15093.0, 15093.0, 15088.0]
    plt.plot(third, label="Solanum Lycopersicum (bad)")

    plt.xlabel("Iteration number")
    plt.ylabel("consensus quality")
    plt.legend()
    plt.savefig("evolutionCOnsensusALL.png")
    plt.close()


def addlabels(x, y, offset):
    for i in range(x):
        plt.text(i, y[i] + offset, int(round(y[i], 0)), rotation=90)


def visualiseExecutionTime(outputDir, timeStart, timeAfterGettingReads,
                           timeAfterConsensus, timeAfterQualityCheck, timeEnd):

    # total time per iteration (also showing identification)
    timeTotal = []
    for i in range(len(timeStart)):
        timeTotal.append(timeEnd[i] - timeStart[i])

    # time quality CHeck per iteration
    timeQuality = []
    for i in range(len(timeStart)):
        timeQuality.append(timeAfterQualityCheck[i] - timeStart[i])

    # time consensus per iteration
    timeConsensus = []
    for i in range(len(timeStart)):
        timeConsensus.append(timeAfterConsensus[i] - timeStart[i])

    # time pre-processing per iteration
    timePreprocess = []
    for i in range(len(timeStart)):
        timePreprocess.append(timeAfterGettingReads[i] - timeStart[i])

    plt.figure()
    plt.title(f"Execution time of each step over time:")

    plt.bar(list(range(len(timeStart))), timeTotal,
            color="y", label="time identification")
    plt.bar(list(range(len(timeStart))), timeQuality,
            color="r", label="time quality check")
    plt.bar(list(range(len(timeStart))), timeConsensus,
            color="g", label="time consensus")
    plt.bar(list(range(len(timeStart))), timePreprocess,
            color="b", label="time pre-processing")

    addlabels(len(timeStart), timeTotal, 0)
    addlabels(len(timeStart), timeQuality, 0)
    addlabels(len(timeStart), timeConsensus, -1)
    addlabels(len(timeStart), timePreprocess, 0.5)

    plt.legend(loc="lower left")

    plt.xlabel("Iteration")
    plt.ylabel("time (s)")
    plt.savefig(f"{outputDir}executionTime.png")
    plt.close()


def saveQualityBothOverTime(outputDir, allqualityreads, allqualityconsensus):

    # now draw quality of consensus over time as quality of read improves
    fig, ax1 = plt.subplots()
    fig.suptitle(
        "Correlation between average read quality and consensus quality:")
    color = 'tab:blue'
    ax1.set_xlabel('iteration number')
    ax1.set_ylabel('average quality of reads in sample', color=color)
    ax1.plot(allqualityreads, color=color)
    ax1.tick_params(axis='y', labelcolor=color)

    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

    # we already handled the x-label with ax1
    color = 'tab:orange'
    ax2.set_ylabel('quality of consensus', color=color)
    ax2.plot(allqualityconsensus, color=color)
    ax2.tick_params(axis='y', labelcolor=color)
    fig.tight_layout()  # otherwise the right y-label is slightly clipped

    plt.savefig(f"{outputDir}evolutionQuality.png")
    plt.close()


def main():

    # visualiseAvgReadQualityTogether()
    # visualiseCOnsQualityTogether()

    pathFastqFile = sys.argv[1]
    geneName = sys.argv[2]

    readLengthDistribution(pathFastqFile, geneName)
    readQualityDistribution(pathFastqFile)


if __name__ == "__main__":
    main()
