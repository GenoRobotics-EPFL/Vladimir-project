from spoa import poa
from src.dataRetriever import *
import os

reads = getReads(11)

# with open("output.fasta", "r") as input:
#     lines = input.readlines()

#     reads = [line[:-1] for line in lines]

consensus, msa = poa(reads, algorithm=0, n=-1, g=-127, e=-127)

# write msa to file
with open("output.fasta", "w+") as result:
    result.writelines(os.linesep.join(msa))

# * create consensus manually

def chooseBestBase(i):

    occurences = {'A': 0, 'C': 0, 'T': 0, 'G': 0}

    for read in msa:

        base = read[i]

        if base == '-':
            continue

        occurences[base] += 1


    # choose base with highest occurence
    highestOccurenceBase = ''
    highestOccurence = 0

    for k, v in occurences.items():
        if v > highestOccurence:
            highestOccurence = v
            highestOccurenceBase = k
        
    print(occurences)

    return highestOccurenceBase





lengthCons = len(msa[0])

cons = ""

for i in range(lengthCons):

    baseChosen = chooseBestBase(i)

    cons += baseChosen

#with open("output.fasta", "w+") as result:
    #result.writelines(os.linesep.join([consensus, cons]))




