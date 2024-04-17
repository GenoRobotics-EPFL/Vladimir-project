import pysam


def getReadsFromFile(pathFastq):

    with pysam.FastxFile(pathFastq, "r") as f:

        sequences = list(f)

        return sequences


def writeReadsToFile(pathFastq, reads):

    with open(pathFastq, "w+") as f:

        readsAsOneString = "\n".join([str(entry) for entry in reads])

        f.write(readsAsOneString)


def writeConsensus(pathCons, cons):

    with open(pathCons, "w+") as f:

        f.write(">Consensus")
        f.write("\n")

        f.write(cons)
        f.write("\n")


def readConsensus(pathCons):

    with open(pathCons, "r") as f:

        lines = f.readlines()

        cons = lines[1]

        cons.replace("\n", "")

        return cons


def writeMsaToFile(pathFile, msa):

    with open(pathFile, "w+") as f:
        f.write("\n".join(msa))
