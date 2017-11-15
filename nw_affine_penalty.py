import numpy as np

delete_matrix = []
insert_matrix = []

def __fill_matrix(chain_a, chain_b, similarity_func, opening_penalty, extension_penalty):
    global delete_matrix, insert_matrix
    matrix = np.zeros((len(chain_a) + 1, len(chain_b) + 1))

    for i in range(1, len(chain_a) + 1):
        matrix[i][0] = matrix[i-1][0] - (opening_penalty + insert_matrix[i-1][0] * extension_penalty)
        insert_matrix[i][0] = insert_matrix[i-1][0] + 1

    for j in range(1, len(chain_b) + 1):
        matrix[0][j] = matrix[0][j-1] - (opening_penalty + delete_matrix[0][j-1] * extension_penalty)
        delete_matrix[i][j] = delete_matrix[i][j-1] + 1

    for i in range(0, len(chain_a)):
        for j in range(0, len(chain_b)):
            match = matrix[i][j] + similarity_func(chain_a[i], chain_b[j])
            delete = matrix[i+1][j] - (opening_penalty + delete_matrix[i+1][j] * extension_penalty)
            insert = matrix[i+1][j] - (opening_penalty + insert_matrix[i][j+1] * extension_penalty)
            matrix[i+1][j+1] = max(match, delete, insert)
            if matrix[i+1][j+1] == delete:
                delete_matrix[i+1][j+1] = delete_matrix[i+1][j] + 1
            elif matrix[i+1][j+1] == insert:
                insert_matrix[i+1][j+1] = insert_matrix[i][j+1] + 1


    return matrix


def __trace_back(chain_a, chain_b, matrix, similarity_func, opening_penalty, extension_penalty):
    alignment_a = ""
    alignment_b = ""
    i = len(chain_a) - 1
    j = len(chain_b) - 1

    while i >= 0 or j >= 0:
        if i >= 0 and j >= 0 and matrix[i + 1][j + 1] == matrix[i][j] + similarity_func(chain_a[i], chain_b[j]):
            alignment_a = chain_a[i] + alignment_a
            alignment_b = chain_b[j] + alignment_b
            i -= 1
            j -= 1
        elif i >= 0 and matrix[i+1][j+1] == matrix[i][j + 1] - (opening_penalty + insert_matrix[i-1][0] * extension_penalty):
            alignment_a = chain_a[i] + alignment_a
            alignment_b = "-" + alignment_b
            i -= 1

        elif j >= 0 and matrix[i+1][j+1] == matrix[i+1][j] - (opening_penalty + delete_matrix[i+1][j] * extension_penalty):
            alignment_a = "-" + alignment_a
            alignment_b = chain_b[j] + alignment_b
            j -= 1

    return alignment_a, alignment_b


def nw_affine_penalty(chain_a, chain_b, similarity_func, opening_penalty, extension_penalty):
    global delete_matrix, insert_matrix
    delete_matrix = np.zeros((len(chain_a) + 1, len(chain_b) + 1))
    insert_matrix = np.zeros((len(chain_a) + 1, len(chain_b) + 1))

    matrix = __fill_matrix(chain_a, chain_b, similarity_func, opening_penalty, extension_penalty)

    aligments = __trace_back(chain_a, chain_b, matrix, similarity_func, opening_penalty, extension_penalty)
    return aligments
