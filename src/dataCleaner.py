"""
This file contains all the preprocessing steps.
The main entry point is "getCleanReads" that is called by the pipelines.

The script can also be directly executed to clean reads from a fastq file and store them back to file:
python3 ./src/dataCleaner.py /pathToFastqFile
"""

import parasail
import sys

import utils
from dataRetriever import *


def removeShortReads(reads, minReadLength):
    return [read for read in reads if len(read.sequence) >= minReadLength]


def removeTooLongReads(reads, maxReadLength):
    return [read for read in reads if len(read.sequence) <= maxReadLength]


def getReferenceReadForOrientation(reads):
    return max(reads, key=lambda read: len(read.sequence))


def findCorrectOrientation(reads, referenceRead=None):
    """ 
    reads now can be the in reverse, and complementary to each other

    from epi2me code wf-amplicon, it is either in correct direction or in reverse AND complementary. Can't be just
    reverse or just complementary.
    """

    def align(s1, s2):
        return parasail.sw_trace_striped_16(s1, s2, open=8, extend=4, matrix=parasail.dnafull)

    if referenceRead is None:
        referenceRead = getReferenceReadForOrientation(reads)

    reference_seq = referenceRead.sequence
    result = [referenceRead]

    for entry in reads:

        if entry == referenceRead:
            continue

        seq = entry.sequence
        fwd_score = align(seq, reference_seq).score

        r = utils.reverse(seq)
        rev_score = align(r, reference_seq).score

        c = utils.complementary(seq)
        comp_score = align(c, reference_seq).score

        rc = utils.complementary(utils.reverse(seq))
        rev_comp_score = align(rc, reference_seq).score

        bestScore = max(fwd_score, rev_comp_score, rev_score, comp_score)

        if fwd_score == bestScore:
            result.append(entry)

        elif rev_score == bestScore:
            entry.sequence = r
            entry.quality = utils.reverse(entry.quality)
            result.append(entry)

        elif comp_score == bestScore:
            entry.sequence = c
            result.append(entry)

        else:
            entry.sequence = rc
            entry.quality = utils.reverse(entry.quality)
            result.append(entry)

    return result, referenceRead


def getQualityRead(read):

    qscores = [ord(ch) - 33 for ch in read.quality]

    meanQScore = sum(qscores) / len(qscores)
    # modeQScore = max(set(qscores), key=qscores.count)

    return meanQScore


def removeBadQualityReads(reads):

    minReadQscore = 10

    goodQualReads = list(
        filter(lambda read: getQualityRead(read) >= minReadQscore, reads))

    return goodQualReads


def removeCorruptedReads(reads):
    """
    When visualising the distribution of base quality, I noticed that somtimes
    bases had quality that were way too high. Those are corrupted reads
    """

    readsCleaned = []

    for read in reads:

        qscores = [ord(ch) - 33 for ch in read.quality]

        if max(qscores) < 80:
            readsCleaned.append(read)

    return readsCleaned


def trimReads(reads):
    """
    https://gigabaseorgigabyte.wordpress.com/2017/06/05/trimming-and-filtering-oxford-nanopore-sequencing-reads/
    """

    # todo: check this perhaps this to too much
    numBasesToTrim = 50

    def trimRead(read):

        read.sequence = read.sequence[numBasesToTrim: -numBasesToTrim]
        read.quality = read.quality[numBasesToTrim: -numBasesToTrim]

        return read

    reads = list(map(trimRead, reads))

    return reads


def getCleanReads(reads, geneName, referenceRead=None):

    print(
        f"Starting cleaning of reads. Original number of reads: \t{len(reads)}")

    minReadLength = 0

    if geneName == "matK":
        avg = 700
        maxReadLength = 900

    elif geneName == "rbcL":
        avg = 750
        maxReadLength = 850

    elif geneName == "psbA-trnH":
        avg = 500
        maxReadLength = 700

    elif geneName == "ITS":
        avg = 700
        maxReadLength = 800
    reads = removeShortReads(reads, minReadLength)
    reads = removeTooLongReads(reads, maxReadLength)

    reads = removeBadQualityReads(reads)

    reads = removeCorruptedReads(reads)

    # EPI2ME also removes the portion of longest reads as those are usually corrupted
    # here I don't do it as sometimes there are too few reads to clean

    # should perhaps also do trimming... but it requires a more complex method than what I did
    # reads = trimReads(reads)

    reads, referenceReadForOrientation = findCorrectOrientation(
        reads, referenceRead)

    print(f"Number of reads after cleanup: \t{len(reads)}")

    return reads, referenceReadForOrientation


def main():

    # do the preprocessing on the original file
    pathInput = sys.argv[1]
    geneName = sys.argv[2]

    reads = getReadsFromFile(pathInput)

    reads, _ = getCleanReads(reads, geneName)

    pathOutput = pathInput + "-CLEANED.fastq"
    writeReadsToFile(pathOutput, reads)


if __name__ == "__main__":
    main()
