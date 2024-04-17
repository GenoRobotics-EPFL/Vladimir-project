import sys
import os
import spoa
import shutil

from utils import *
from dataRetriever import *


def _getDraftConsensus(pathToReads):
    print("making draft consensus using spoa")

    reads = getReadsFromFile(pathToReads)
    sequences = [read.sequence for read in reads]

    min_cov = int(round(len(sequences) * 0.15))

    # run spoa will all sequence first
    cons, _ = spoa.poa(sequences, min_coverage=min_cov, genmsa=False)

    # run spoa a second time with the previous result as first read
    cons, _ = spoa.poa([cons, *sequences], min_coverage=min_cov, genmsa=False)

    # write consensus to file
    pathOutput = os.path.join("outputSpoa", "draftConsensusSpoa.fasta")
    writeConsensus(pathOutput, cons)

    return cons


def _getConsensusMedaka(pathToReads):
    print("making final consensus usind medaka")

    pathToDraft = os.path.join("outputSpoa", "draftConsensusSpoa.fasta")
    pathToOutputDir = os.path.join("outputMedaka")

    medakaCommand = f"medaka_consensus -i {pathToReads} -d {pathToDraft} -o {pathToOutputDir} -m r941_min_high_g303"
    os.system(medakaCommand)

    # read medaka output
    pathToOutputConsensus = os.path.join(pathToOutputDir, "consensus.fasta")
    cons = readConsensus(pathToOutputConsensus)

    return cons


def getConsensus(pathToReads):
    """
    get a conensus of multiple sequences (bases only), using the spoa library.
    the consensus will compare all sequences together to find the best match.
    """

    # must clean directory in case it exists. Medaka complains otherwise
    pathToOutputSpoa = os.path.join("outputSpoa")
    pathToOutputMedaka = os.path.join("outputMedaka")

    dirToClear = os.path.join(pathToOutputSpoa, "*")
    os.system(f"rm {dirToClear}")
    dirToClear = os.path.join(pathToOutputMedaka, "*")
    os.system(f"rm {dirToClear}")

    # first get the consensus of spoa as draft consensus
    _getDraftConsensus(pathToReads)

    # use medaka to improve result, by remapping the reads
    cons = _getConsensusMedaka(pathToReads)

    return cons


def main():

    pathToReads = sys.argv[1]

    cons = getConsensus(pathToReads)

    print(f"\nFinal cons:\n\n{cons} \n")


if __name__ == "__main__":
    main()
