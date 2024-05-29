"""
This file takes as input a whole fastq file, and will simulate
the output of the sequencing as if it is happening on the computer.

the calculations are based on report_FAV14934_20221125_1248_ee7626ed.html

usage:

python simulateRealTimeData.py path/to/fastqfile
"""

import sys
import os
import shutil
import time

# total number of passing reads after 4 hours of sequencing
total_num_reads_atEnd = 200_000

# number of 10 minute intervals in the time all the reads were generated
num_intervals = 4 * 6  # 4 hours, each with 6 intervals

# number of passing reads received each 10 minutes
num_reads_each_iteration = total_num_reads_atEnd / num_intervals

# number of genes/barcodes all these reads have to go towards
num_genes = 48

# number of reads per interation that are for a given gene
num_reads_per_gene_per_iteration = int(num_reads_each_iteration / num_genes) + 1

outputDir = "fastqpass"


def simulateOutput(pathFastq: str):

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
            nameFile = "iteration" + str(iterationNum) + ".fastq"
            outputFilePath = os.path.join(outputDir, nameFile)

            with open(outputFilePath, "w+") as outputFile:
                outputFile.writelines(currIterationContent)

            # time.sleep(1 * 60 * 10)  # sleep 10 minutes
            iterationNum += 1


def main():

    pathFastq = sys.argv[1]

    # create directory where output is stored
    if os.path.exists(outputDir):
        shutil.rmtree(outputDir)

    os.makedirs(outputDir)

    # start simulation
    simulateOutput(pathFastq)


if __name__ == "__main__":
    main()
