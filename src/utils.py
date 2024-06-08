"""
Util script with methods that other files can use
"""

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