
def complementary(seq):
    """
    get the complementary of a sequence of bases
    """

    complementaryBase = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

    return "".join([complementaryBase[base] for base in seq])

def reverse(seq):
    """
    get the reverse of the given list
    """

    return seq[::-1]

def writeStringToFile(pathFile, str):

    with open(pathFile, "w+") as f:
        f.write(str)