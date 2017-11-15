import numpy as np

from blosum62 import blosum_62_scoring
from smith_waterman import smith_waterman
from hirschberg import hirschberg
from needleman_wunsch import needleman_wunsch
from nw_affine_penalty import nw_affine_penalty


def main():
    f = open('test_input_1.txt', 'r')
    str = f.read()
    chain_a_raw, chain_b_raw = str.split("\n")
    chain_a = np.array(list(chain_a_raw))
    chain_b = np.array(list(chain_b_raw))
    # alignments = needleman_wunsch(chain_a, chain_b, blosum_62_scoring, gap_penalty=-5)
    # alignments = hirschberg(chain_a, chain_b, -5, blosum_62_scoring)
    # alignments = smith_waterman(chain_a, chain_b, blosum_62_scoring, gap_penalty=-5)
    alignments = nw_affine_penalty(chain_a, chain_b, blosum_62_scoring, 5, 2)
    print("Result\n", alignments[0], "\n", alignments[1])
    f.close()


if __name__ == "__main__":
    main()
