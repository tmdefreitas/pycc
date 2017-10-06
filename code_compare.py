# code_compare.py

from tokenize import tokenize
import token
from io import BytesIO
import numpy as np


def _tokenize(s):
    for tok in tokenize(BytesIO(s.encode('utf-8')).readline):
        yield tok


def align(s1, s2, sim_func, gap_penalty=-1, method="global"):
    m = np.zeros((len(s1)+1, len(s2)+1))


    _max_score = 0
    _max_i = 0
    _max_j = 0

    # initialize first row and column, if doing global
    if method == "global":
        for i in range(len(s1)):
            m[i+1][0] = m[i][0] + gap_penalty

        for j in range(len(s2)):
            m[0][j+1] = m[0][j] + gap_penalty

    # fill in the table
    for i in range(1,len(s1)+1):
        for j in range(1, len(s2)+1):
            # m[i][j] is the max of 

            left = m[i][j-1] + gap_penalty
            up   = m[i-1][j] + gap_penalty
            diag = m[i-1][j-1] + sim_func(s1[i-1], s2[j-1]) # compare current characters, remembering i and j index the table, one less than the strings

            score = max(left, up, diag)

            if method == "local" and score < 0:
                score = 0

            # keep track of the maximum score position, for local alignments
            if score > _max_score:
                _max_i = i
                _max_j = j
                _max_score = score

            m[i][j] = score

    # Backtrack to find the 
    aln_s1 = []
    aln_s2 = []

    if method == "global":
        b_i = len(s1)
        b_j = len(s2)
    else:
        b_i = _max_i
        b_j = _max_j

    while b_i > 0 and b_j > 0:
        # check to see which box this came from
        left = m[b_i][b_j-1] + gap_penalty
        up   = m[b_i-1][b_j] + gap_penalty
        diag = m[b_i-1][b_j-1] + sim_func(s1[b_i-1], s2[b_j-1]) 

        # find out which sequnce to take elements from. None indicates a gap
        if m[b_i][b_j] == diag:
            aln_s1.insert(0, s1[b_i-1])
            aln_s2.insert(0, s2[b_j-1])
            b_i -= 1
            b_j -= 1
        elif m[b_i][b_j] == left:
            aln_s1.insert(0, None)
            aln_s2.insert(0, s2[b_j-1])
            b_j -= 1
        else:
            aln_s1.insert(0, s1[b_i-1])
            aln_s2.insert(0, None)
            b_i -= 1

        # end local alignments early
        if method == "local" and m[b_i][b_j] == 0:
            break

    ret_score = _max_score if method=="local" else m[-1][-1]

    return ret_score, aln_s1, aln_s2, m


def simple_sim(t1, t2):
    return 1 if t1 == t2 else -3

# Prints an alignment of sequences, where seq2str is a function that turns an element of seq into a printable string
def print_seq_alignment(seq1, seq2, seq2str, max_line_length=80):
    assert(len(seq1) == len(seq2))

    outline1 = ""
    outline2 = ""

    strings1 = [seq2str(c) for c in seq1]
    strings2 = [seq2str(c) for c in seq2]

    for i in range(len(seq1)):
        s1 = strings1[i]
        s2 = strings2[i]

        if len(s1) > len(s2):
            s2 = s2.ljust(len(s1))
        else:
            s1 = s1.ljust(len(s2))

        outline1 += s1
        outline2 += s2

        if len(outline1) > max_line_length:
            print(outline1)
            print(outline2 + "\n")
            outline1 = ""
            outline2 = ""

    if outline1:
        print(outline1)
        print(outline2 + "\n")


def charseq2str(c):
    if c:
        return c
    else:
        return "-"

sentence1 = "This is a long sentence with words in it."
sentence2 = "This is a shorter sentence."

aln = align(sentence1, sentence2, simple_sim)

print(aln[0])
print_seq_alignment(aln[1], aln[2], charseq2str)



