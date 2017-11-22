import numpy as np


def __fill_matrix(chain_a, chain_b, similarity_func, gap_penalty = -1):
    matrix = np.zeros((len(chain_a) + 1, len(chain_b) + 1))

    for i in range(len(chain_a) + 1):
        matrix[i][0] = gap_penalty * i

    for j in range(len(chain_b) + 1):
        matrix[0][j] = gap_penalty * j

    for i in range(0, len(chain_a)):
        for j in range(0, len(chain_b)):
            match = matrix[i][j] + similarity_func(chain_a[i], chain_b[j])
            delete = matrix[i][j + 1] + gap_penalty
            insert = matrix[i + 1][j] + gap_penalty
            matrix[i + 1][j + 1] = max(match, delete, insert)

    return matrix


def __trace_back(chain_a, chain_b, matrix, similarity, gap_penalty = -1):
    alignment_a = ""
    alignment_b = ""
    i = len(chain_a) - 1
    j = len(chain_b) - 1
    while i >= 0 or j >= 0:
        if i >= 0 and j >= 0 and matrix[i + 1][j + 1] == matrix[i][j] + similarity(chain_a[i], chain_b[j]):
            alignment_a = chain_a[i] + alignment_a
            alignment_b = chain_b[j] + alignment_b
            i -= 1
            j -= 1
        elif i >= 0 and matrix[i+1][j+1] == matrix[i][j + 1] + gap_penalty:
            alignment_a = chain_a[i] + alignment_a
            alignment_b = "-" + alignment_b
            i -= 1

        else:
            alignment_a = "-" + alignment_a
            alignment_b = chain_b[j] + alignment_b
            j -= 1

    return alignment_a, alignment_b


def needleman_wunsch(chain_a, chain_b, similarity_func, gap_penalty = -1, print_matrix=False):
    matrix = __fill_matrix(chain_a, chain_b, similarity_func, gap_penalty=gap_penalty)
    if print_matrix:
        print(len(matrix[0]))

    aligments = __trace_back(chain_a, chain_b, matrix, similarity_func, gap_penalty=gap_penalty)
    score = matrix[-1][-1]
    return aligments, score
