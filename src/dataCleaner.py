import parasail
import sys

import utils
from dataRetriever import *


def removeShortReads(reads):

    minReadLength = 300

    return [read for read in reads if len(read.sequence) >= minReadLength]


def removeTooLongReads(reads, maxReadLength):
    return [read for read in reads if len(read.sequence) <= maxReadLength]


def getReferenceReadForOrientation(reads):
    return max(reads, key=lambda read: len(read.sequence))


def findCorrectOrientation(reads, referenceRead=None):
    """ 
    reads now can be the in reverse, and complementary to each other

    from epi2me code wf-amplicon, it is iehter in correct direction or in reverse AND complementary. Can't be just
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
    modeQScore = max(set(qscores), key=qscores.count)

    return meanQScore


def removeBadQualityReads(reads, minReadQscore):

    goodQualReads = list(
        filter(lambda read: getQualityRead(read) >= minReadQscore, reads))

    return goodQualReads


def trimReads(reads):
    """
    https://gigabaseorgigabyte.wordpress.com/2017/06/05/trimming-and-filtering-oxford-nanopore-sequencing-reads/
    """

    numBasesToTrim = 50

    def trimRead(read):

        read.sequence = read.sequence[numBasesToTrim: -numBasesToTrim]
        read.quality = read.quality[numBasesToTrim: -numBasesToTrim]

        return read

    reads = list(map(trimRead, reads))

    return reads


def getCleanReads(reads, maxReadLength, minReadQscore,  referenceRead=None):

    print(f"Original number of reads: \t{len(reads)}")

    reads = removeShortReads(reads)
    reads = removeTooLongReads(reads, maxReadLength)

    reads = removeBadQualityReads(reads, minReadQscore)

    # EPI2ME also removes the portion of longest reads as those are usually corrupted
    # here I don't do it as sometimes there are too few reads to clean

    reads = trimReads(reads)

    reads, referenceReadForOrientation = findCorrectOrientation(reads, referenceRead)

    print(f"Number of reads after cleanup: \t{len(reads)}")

    return reads, referenceReadForOrientation 




def main():

    # do the preprocessing on the original file
    pathInput = sys.argv[1]
    pathOutput = pathInput + str("_cleaned.fastq")

    reads = getReadsFromFile(pathInput)

    maxReadLength = 750
    minReadQscore = 10
    reads, _ = getCleanReads(reads, maxReadLength, minReadQscore)

    writeReadsToFile(pathOutput, reads)


if __name__ == "__main__":
    main()
