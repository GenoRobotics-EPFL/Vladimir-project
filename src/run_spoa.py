"""Interleave forward and reverse reads to improve SPOA performance."""

import sys
import os

import parasail
import pysam
import spoa

from preProcess import *
from dataRetriever import *


def align(s1, s2):
    """Use parasail to align two sequences.

    :param s1: string with first DNA sequence
    :param s2: string with second DNA sequence
    :return: parasail alignment Result object
    """

    return parasail.sw_trace_striped_16(
        s1, s2, open=8, extend=4, matrix=parasail.dnafull
    )


def interleave_lists(l1, l2):
    """Uniformly interleave two lists.

    For example:
    ```
    >>> interleave_lists([1, 2, 3, 4], list('ABCDEF'))
    [1, 'A', 'B', 2, 'C', 3, 'D', 'E', 4, 'F']
    ```

    :param l1: first list
    :param l2: second list
    :return: list containing items of `l1` and `l2` uniformly interleaved
    """
    interleaved = []
    rate_1 = 1 / len(l1)
    rate_2 = 1 / len(l2)

    c_1, c_2 = 0, 0

    itr_1 = iter(l1)
    itr_2 = iter(l2)

    while True:
        try:
            if c_1 > c_2:
                interleaved.append(next(itr_2))
                c_2 += rate_2
            else:
                interleaved.append(next(itr_1))
                c_1 += rate_1
        except StopIteration:
            break

    return interleaved


def main():

    pathfastq = os.path.join(sys.argv[1])

    reads, _ = getCleanReads(pathfastq)

    # read input file and determine the orientation of the reads
    fwd = []
    rev = []


    # this is our reference read
    # todo: take the longest read instead as reference !
    first_seq = reads[0].sequence
    fwd.append(first_seq)

    for entry in reads[1:]:

        seq = entry.sequence
        rc = complementary(reverse(seq))

        fwd_score = align(seq, first_seq).score
        rev_score = align(rc, first_seq).score

        if fwd_score > rev_score:
            fwd.append(seq)
        else:
            rev.append(rc)

    # uniformly interleave the reads
    interleaved_reads = interleave_lists(fwd, rev)


    # run spoa
    min_cov = int(round(len(interleaved_reads) * 0.15))
    cons, _ = spoa.poa(interleaved_reads, min_coverage=min_cov, genmsa=False)

    # run spoa a second time with the previous result as first read
    cons, _ = spoa.poa([cons, *interleaved_reads], min_coverage=min_cov, genmsa=False)

    # write out the result
    # sys.stdout.write(os.linesep.join(msa))

    print(cons)

if __name__ == "__main__":
    main()