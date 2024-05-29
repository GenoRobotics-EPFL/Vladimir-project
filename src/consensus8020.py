import sys
import os
import spoa

from utils import *
from dataRetriever import *

def getDraftConsensus(pathToReads):
    """
    get a conensus of multiple sequences (bases only), using the spoa library.
    the consensus will compare all sequences together to find the best match.
    """

    print("making draft consensus using spoa")

    reads = getReadsFromFile(pathToReads)
    sequences = [read.sequence for read in reads]

    min_cov = int(round(len(sequences) * 0.15))

    # run spoa will all sequence first
    cons, _ = spoa.poa(sequences, min_coverage=min_cov, genmsa=False)

    # run spoa a second time with the previous result as first read
    cons, _ = spoa.poa([cons, *sequences], min_coverage=min_cov, genmsa=False)

    # write consensus to file
    pathOutput = os.path.join("outputSpoa", "consensusSpoa.fasta")
    writeConsensus(pathOutput, cons)

    return cons

def getConsensusMedaka(pathToReads, pathToDraft):
    print("making final consensus usind medaka")

    pathToOutputDir = os.path.join("outputMedaka")
    medakaCommand = f"medaka_consensus -i {pathToReads} -d {pathToDraft} -o {pathToOutputDir}"

    medakaCommand = f"medaka_consensus -q -i {pathToReads} -d {pathToDraft} -o {pathToOutputDir} -m r941_min_high_g303 > /dev/null 2> /dev/null"
    os.system(medakaCommand)

    # read medaka output
    pathToOutputConsensus = os.path.join(pathToOutputDir, "consensus.fastq")
    cons = getReadsFromFile(pathToOutputConsensus)[0]

    return cons

def getConsensusAndQuality(pathToReads):

    # must clean directory in case it exists. Medaka complains otherwise
    os.system("rm ./outputSpoa/*")
    os.system("rm ./outputMedaka/*")

    # first get the consensus of spoa as draft consensus
    getDraftConsensus(pathToReads)

    # use medaka to improve result, by remapping the reads
    pathToDraft = os.path.join("outputSpoa", "consensusSpoa.fasta")
    cons = getConsensusMedaka(pathToReads, pathToDraft)

    return cons

def getConsensus(pathToReads):

    cons = getConsensusAndQuality(pathToReads)

    return cons.sequence

def splitReads(reads):

    reads.sort(key=lambda read: len(read.sequence), reverse=True)

    indexTwenty = int(len(reads) * 0.2)

    return reads[:indexTwenty], reads[indexTwenty:]


def main():

    pathToReads = sys.argv[1]
    reads = getReadsFromFile(pathToReads)

    twenty, eighty = splitReads(reads)

    pathToReadsTwenty = pathToReads + "-twenty"
    writeReadsToFile(pathToReadsTwenty, twenty)
    pathToReadsEighty = pathToReads + "-eighty"
    writeReadsToFile(pathToReadsEighty, eighty)

    cons20 = getConsensus(pathToReadsTwenty)

    pathToDraft = pathToReads + "-cons20"
    writeConsensus(pathToDraft, cons20)

    cons80 = getConsensusMedaka(pathToReadsEighty, pathToDraft)

    print(f"\nFinal cons:\n\n{str(cons80)} \n")


if __name__ == "__main__":
    main()