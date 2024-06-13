"""
This file is used to retrieve and write reads from/to fastq files.
"""

import pysam
import subprocess
import time
import os


def getReadsFromFile(pathFastq):

    with pysam.FastxFile(pathFastq, "r") as f:

        sequences = list(f)

        return sequences


def waitForReadsFromFile(pathFastq):

    # wait for the file to come to existance
    while not os.path.exists(pathFastq):
        print(f"Waiting for file: {pathFastq}")
        time.sleep(10)

    with pysam.FastxFile(pathFastq, "r") as f:

        sequences = list(f)

        return sequences


def waitForNewReadsFile(readDir, knownFileNames):

    while True:

        output = subprocess.check_output(["ls", readDir])
        fileNames = output.decode('utf-8').split("\n")

        newReads = []
        for fileName in fileNames:
            if fileName == "":
                continue

            if fileName not in knownFileNames:
                print(f"Taking new reads in file: {fileName}")

                knownFileNames.add(fileName)

                reads = getReadsFromFile(readDir + fileName)
                newReads += reads

        if len(newReads) == 0:
            print("Waiting for new reads")
            time.sleep(10)
            continue

        return newReads, knownFileNames


def writeReadsToFile(pathFastq, reads):

    with open(pathFastq, "w+") as f:

        readsAsOneString = "\n".join([str(entry) for entry in reads])

        f.write(readsAsOneString + "\n")


def writeConsensus(pathCons, cons):

    with open(pathCons, "w+") as f:

        f.write(">Consensus")
        f.write("\n")

        f.write(cons)
        f.write("\n")


def readConsensus(pathCons):

    with open(pathCons, "r") as f:

        lines = f.readlines()

        cons = lines[1]

        cons.replace("\n", "")

        return cons


def writeMsaToFile(pathFile, msa):

    with open(pathFile, "w+") as f:
        f.write("\n".join(msa))


def getLinesFromFile(pathFile):

    with open(pathFile, "r") as f:

        lines = f.readlines()

        lines = list(filter(lambda line: line != "", lines))

        lines = list(map(lambda line: line.replace("\n", ""), lines))

        return lines
