import numpy as np


def __fill_matrix(chain_a, chain_b, similarity_func, gap_penalty = -1):
    matrix = np.zeros((len(chain_a) + 1, len(chain_b) + 1))

    for i in range(0, len(chain_a)):
        for j in range(0, len(chain_b)):
            match = matrix[i][j] + similarity_func(chain_a[i], chain_b[j])
            delete = matrix[i][j + 1] + gap_penalty
            insert = matrix[i + 1][j] + gap_penalty
            matrix[i + 1][j + 1] = max(match, delete, insert, 0)

    return matrix


def __maximum_score(matrix):
    n, m = matrix.shape

    max_val = 0
    max_val_i = n - 1
    max_val_j = m - 1

    for i in range(n-1, 0, -1):
        for j in range(m-1, 0, -1):
            if matrix[i][j] > max_val:
                max_val = matrix[i][j]
                max_val_i = i
                max_val_j = j

    return max_val, max_val_i, max_val_j

def __trace_back(chain_a, chain_b, matrix, similarity, gap_penalty = -1):

    max_val, max_val_i, max_val_j = __maximum_score(matrix)

    alignment_a = ""
    alignment_b = ""
    i = max_val_i
    j = max_val_j
    value = max_val
    while (i > 0 or j > 0) and value != 0:
        if i > 0 and j > 0 and matrix[i][j] == matrix[i-1][j-1] + similarity(chain_a[i - 1], chain_b[j - 1]):
            alignment_a = chain_a[i - 1] + alignment_a
            alignment_b = chain_b[j - 1] + alignment_b
            i -= 1
            j -= 1
        elif i > 0 and matrix[i][j] == matrix[i - 1][j] + gap_penalty:
            alignment_a = chain_a[i] + alignment_a
            alignment_b = "-" + alignment_b
            i -= 1

        else:
            alignment_a = "-" + alignment_a
            alignment_b = chain_b[j - 1] + alignment_b
            j -= 1

        value = matrix[i][j]

    return alignment_a, alignment_b

def smith_waterman(chain_a, chain_b, similarity_func, gap_penalty = -1, print_matrix=False):
    matrix = __fill_matrix(chain_a, chain_b, gap_penalty=gap_penalty, similarity_func=similarity_func)
    if print_matrix:
        print(len(matrix[0]))

    alignments = __trace_back(chain_a, chain_b, matrix, gap_penalty=gap_penalty, similarity=similarity_func)
    score = matrix[-1][-1]
    return alignments, score