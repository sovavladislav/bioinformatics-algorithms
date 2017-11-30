from needleman_wunsch import needleman_wunsch
from blosum62 import blosum_62_scoring
import numpy as np

def msa(chains, gap_penalty):
    matrix = np.zeros((len(chains), len(chains)))
    #compute pairwise alignment distance between every pair
    for i in range(len(chains)):
        for j in range(i+1, len(chains)):
            _, matrix[i][j] = needleman_wunsch(chains[i], chains[j], score, gap_penalty)
            matrix[j][i] = matrix[i][j]

    # choose minimum SUM distance of chains and denote it as center chain
    min_row = matrix.sum(axis=0).argmin()

    # set this chain first
    move_elem(chains, min_row, 0)

    #multiple aligment by pairwise alignment
    for i in range(1, len(chains)-1):
        for j in range(1, i+1):
            alignments, s = needleman_wunsch(chains[0], chains[j], score, gap_penalty)
            chains[0], chains[j] = alignments[0], alignments[1]
        alignments, s = needleman_wunsch(chains[0], chains[i+1], score, gap_penalty)
        chains[0], chains[i+1] = alignments[0], alignments[1]

    for i in range(1, len(chains)):
        alignments, s = needleman_wunsch(chains[0], chains[i], score, gap_penalty)
        chains[0], chains[i] = alignments[0], alignments[1]

    #return back chain
    move_elem(chains, 0, min_row)

    res = []
    for i in chains:
        res.insert(len(res), ''.join(i))

    return res

def score(a,b):
    if a == '-' or b == '-':
        return -5
    else:
        return blosum_62_scoring(a,b)

def move_elem(chains, i, j):
    c = chains[i]
    del chains[i]
    chains.insert(j, c)