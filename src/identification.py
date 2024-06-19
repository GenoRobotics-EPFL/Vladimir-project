"""
This file executes the blastn queries to detect a plant using a consensus
"""

import sys
import os
from Bio.Blast import NCBIXML


class NoIdentification(Exception):
    pass


def identify(pathConsensus, db):

    # clean blastn output file
    outputDir = "./outputBlastn"
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)
    else:
        os.system(f"rm {outputDir}/*")

    pathOutput = "outputBlastn/temp.xml"

    blastnCommand = f"blastn -query {pathConsensus} -db {db} -out {pathOutput} -outfmt 5 -max_hsps 1"
    os.system(blastnCommand)

    # now parse results
    blast_records = NCBIXML.parse(open(pathOutput, "r"))
    blast_record = next(blast_records)  # you only have one query

    if len(blast_record.alignments) == 0:
        raise NoIdentification()

    alignment = blast_record.alignments[0]  # you take the best alignment

    sumcoverage = 0
    sumidentity = 0
    sumAlignLenth = 0

    for hsp in alignment.hsps:
        if hsp.expect < 0.01:
            sumcoverage += (hsp.query_end - hsp.query_start + 1)

            sumidentity += hsp.identities
            sumAlignLenth += hsp.align_length

    coverage = (sumcoverage / blast_record.query_length) * 100
    identity = (sumidentity / sumAlignLenth) * 100

    return (alignment.hit_def, coverage, identity)


def printResultIdentification(result):

    print(
        f"The species name is: \n {result[0]}\ncoverage: {result[1]}\nidentity: {result[2]}")


def main():

    pathConsensus = os.path.join(sys.argv[1])

    db = sys.argv[2]

    result = identify(pathConsensus, db)

    printResultIdentification(result)


if __name__ == "__main__":
    main()
