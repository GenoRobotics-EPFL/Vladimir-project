from dataRetriever import *
from consensusIterative import *
from readpaf import parse_paf


def getReadWithName(reads, queryName):

    for read in reads:
        if read.name == queryName:
            return read


def firstIteration(pathFastq):
    print("first iteration starting")

    reads = getReadsFromFile(pathFastq)
    sequences = [entry.sequence for entry in reads]

    # this will guide for the temporary consensus
    tempCons = spoaCons(sequences)
    tempConsPath = os.path.join("./tempCons.fasta")
    writeConsensus(tempConsPath, tempCons)  # you must do this

    # align the original reads on the consensus with minimap
    pathResultMinimap = os.path.join("./tempMinimap.paf")
    os.system(
        f"minimap2 -x map-ont {tempConsPath} {pathFastq} > {pathResultMinimap}")

    # recreate consensus from minimap result
    tempConsObj = TemporaryConsensus(tempCons)

    # align minimap output to the temporary consensus
    for entry in parse_paf(open(pathResultMinimap, "r")):

        read = getReadWithName(reads, entry.query_name)

        tempConsObj.acceptRead(
            read, entry.query_start, entry.query_end, entry.target_start, entry.target_end)

    print(tempConsObj.getCons())


def laterIteration(iterationNum):
    print("later iteration starting")


def start():

    iterationNum = 0

    while True:

        if iterationNum == 0:
            firstIteration()
        else:
            laterIteration(iterationNum)

        iterationNum += 1


def main():
    print("starting real time streaming")

    pathFastq = sys.argv[1]

    firstIteration(pathFastq)


if __name__ == "__main__":
    main()
