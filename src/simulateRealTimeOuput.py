"""
This file takes as input a whole fastq file, and a time between iterations, and will simulate
the output of the sequencing as if it is happening on the computer.

the calculations are based on report_FAV14934_20221125_1248_ee7626ed.html

usage:

python simulateRealTimeData.py path/to/fastqfile timeBetweenIteration

where timeBetween iterations is in minutes
"""

import sys
import os
import shutil
import time

outputDir = "./fastqpass/"
if not os.path.exists(outputDir):
    os.makedirs(outputDir)
else:
    os.system(f"rm {outputDir}*")

# total number of reads after 4 hours of sequencing
total_num_reads_atEnd = 200_000

# number of time between each iteration where the sequencer created a new fastq files
timeBetweenIter = 1 

# number of intervals in the time all the reads were generated
num_intervals = (4 * 60) / timeBetweenIter

# number of passing reads received each iteration
num_reads_each_iteration = total_num_reads_atEnd / num_intervals

# number of genes/barcodes all these reads have to go towards
num_genes = 48

# number of reads per interation that are for a given gene
num_reads_per_gene_per_iteration = int(
    num_reads_each_iteration / num_genes) + 1

def simulateOutput(pathFastq, sleepBetweenIter):

    with open(pathFastq, "r") as fastq:

        iterationNum = 0

        reachedEOF = False
        while not reachedEOF:

            # read following amount of reads
            # each reach is written on 4 lines
            currIterationContent = []
            for _ in range(num_reads_per_gene_per_iteration * 4):

                line = fastq.readline()

                if line == "":  # check if you reached the end of the file
                    reachedEOF = True
                    break
                else:
                    currIterationContent.append(line)

            # write to new file
            nameFile = f"iteration{iterationNum}.fastq"
            outputFilePath = os.path.join(outputDir, nameFile)

            with open(outputFilePath, "w+") as outputFile:
                outputFile.writelines(currIterationContent)

            time.sleep(sleepBetweenIter * 60)
            iterationNum += 1


def main():

    pathFastq = sys.argv[1]
    timeBetweenIter = int(sys.argv[2])  # in minutes

    # create directory where output is stored
    if os.path.exists(outputDir):
        shutil.rmtree(outputDir)
    os.makedirs(outputDir)

    # start simulation
    simulateOutput(pathFastq, timeBetweenIter)


if __name__ == "__main__":
    main()
