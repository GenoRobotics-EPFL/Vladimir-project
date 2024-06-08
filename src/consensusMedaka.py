"""
THis file is used to get a normal consensus usind Medaka and SPOA.
no preprocessing is done on the reads

Usage:
pyhton3 src/consensusMedaka.py ./allData/fastqFile
"""

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

    print("making consensus using spoa")
    os.system("rm ./outputSpoa/*")

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
    print("making consensus using medaka")

    pathToOutputDir = "./outputMedaka"
    os.system("rm ./outputMedaka/*")

    medakaCommand = f"medaka_consensus -q -i {pathToReads} -d {pathToDraft} -o {pathToOutputDir} -m r941_min_high_g303 > /dev/null 2> /dev/null"
    os.system(medakaCommand)

    # read medaka output
    pathToOutputConsensus = os.path.join(pathToOutputDir, "consensus.fastq")
    outputReads = getReadsFromFile(pathToOutputConsensus)

    if len(outputReads) == 0:  # it means medaka couldn't improve the draft consensus, so you just return the draft consensus
        print("medaka couldn't improve consensus")
        return getReadsFromFile(pathToDraft)[0]

    return outputReads[0]


def getConsensusAndQuality(pathToReads):

    # first get the consensus of spoa as draft consensus
    getDraftConsensus(pathToReads)

    # use medaka to improve result, by remapping the reads
    pathToDraft = os.path.join("outputSpoa", "consensusSpoa.fasta")
    cons = getConsensusMedaka(pathToReads, pathToDraft)

    return cons


def getConsensus(pathToReads):

    cons = getConsensusAndQuality(pathToReads)

    return cons.sequence


def getConsensus8020(best20Reads, allReads):

    # first get the consensus of spoa as draft consensus
    getDraftConsensus(best20Reads)

    # use medaka to improve result, by remapping the reads
    pathToDraft = os.path.join("outputSpoa", "consensusSpoa.fasta")
    cons = getConsensusMedaka(allReads, pathToDraft)

    return cons.sequence


def main():

    pathToReads = sys.argv[1]

    cons = getConsensusAndQuality(pathToReads)

    print(f"\nFinal cons:\n\n{str(cons)} \n")


if __name__ == "__main__":
    main()
