""" 

downgrade a fastq file to make the reads a bit worse

"""

import sys
import random


from dataRetriever import *
from dataCleaner import *
from visualise import *


def removeQualityBase(reads):

    for read in reads:

        newQuality = ""

        for quality in read.quality:
            newQuality += chr(max(ord(quality) - 6, 33))

        read.quality = newQuality

    return reads


def adaptBaseToScore(reads):

    otherBases = {'A': ['C', 'T', 'G'], 'C': ['A', 'G', 'T'],
                  'T': ['A', 'C', 'G'], 'G': ['A', 'C', 'T']}

    for read in reads:

        newSequence = ""
        modifications = 0

        for i, quality in enumerate(read.quality):

            score = ord(quality) - 33
            probError = 10 ** (- score / 10)

            if random.random() < probError:
                newSequence += random.choice(otherBases[read.sequence[i]])
                modifications += 1
            else:
                newSequence += read.sequence[i]

        read.sequence = newSequence

    return reads


def cutReadsInTwo(reads):

    cutReadProb = 0.3

    for read in reads:

        if random.random() < cutReadProb:
            cutIndex = random.randint(0, len(read.sequence))
            read.sequence = read.sequence[:cutIndex]
            read.quality = read.quality[:cutIndex]

    return reads


def downgradeReads(pathToFastq):

    reads = getReadsFromFile(pathToFastq)

    reads = removeQualityBase(reads)
    reads = adaptBaseToScore(reads)

    reads = cutReadsInTwo(reads)

    writeReadsToFile(pathToFastq + "-DOWNGRADED.fastq", reads)


def main():

    pathToFastq = sys.argv[1]

    downgradeReads(pathToFastq)


if __name__ == "__main__":
    main()
