import sys
import os

from utils import *
from dataRetriever import *


def getConsensusMinimap(pathToReads, pathToDraft):
    
    pathToMinimapResult = "./outputMinimap/minimapResult.paf"
    command = f"minimap2 -x map-ont {pathToDraft} {pathToReads} > {pathToMinimapResult}"
    os.system(command)

    pathToRaconResult = "./outputRacon/raconResult.fasta"
    command = f"racon -m 8 -x -6 -g -8 -w 500 {pathToReads} {pathToMinimapResult} {pathToDraft} > {pathToRaconResult}"
    os.system(command)

    consensus = getReadsFromFile(pathToRaconResult)[0].sequence

    return consensus



def main():

    pathToReads = sys.argv[1]
    pathToDraft = sys.argv[2]

    consensus = getConsensusMinimap(pathToReads, pathToDraft)

    print("\n\n" + consensus)

if __name__ == "__main__":
    main()