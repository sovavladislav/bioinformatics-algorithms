import numpy as np

from needleman_wunsch import needleman_wunsch


def __nw_score_last_line(chain_a, chain_b, similarity_func, gap_penalty):
    len_chain_a = len(chain_a)
    len_chain_b = len(chain_b)
    previous_row = np.zeros(shape=(len_chain_b + 1), dtype=np.int)
    current_row = np.zeros(shape=(len_chain_b + 1), dtype=np.int)

    for j in range(1, len_chain_b + 1):
        previous_row[j] = previous_row[j - 1] + gap_penalty

    for i in range(1, len_chain_a + 1):
        current_row[0] = gap_penalty + previous_row[0]
        for j in range(1, len_chain_b + 1):
            score_sub = previous_row[j - 1] + similarity_func(chain_a[i - 1], chain_b[j - 1])
            score_del = previous_row[j] + gap_penalty
            score_ins = current_row[j - 1] + gap_penalty
            current_row[j] = max(score_sub, score_del, score_ins)

        previous_row = current_row
        current_row = [0] * (len_chain_b + 1)

    return previous_row


def __hirschberg(chain_a, chain_b, similarity_func, gap_penalty):
    first_aligned, second_aligned = [], []
    length_chain_a, length_chain_b = len(chain_a), len(chain_b)

    if length_chain_a == 0:
        for i in range(length_chain_b):
            first_aligned.append("-" * len(chain_b[i]))
            second_aligned.append(chain_b[i])
    elif length_chain_b == 0:
        for i in range(length_chain_a):
            first_aligned.append(chain_a[i])
            second_aligned.append("-" * len(chain_a[i]))

    elif length_chain_a == 1 or length_chain_b == 1:
        first_aligned, second_aligned = needleman_wunsch(chain_a, chain_b, similarity_func, gap_penalty)

    else:

        # Divide and Conquer

        middle_chain_a = int(length_chain_a / 2)

        row_left = __nw_score_last_line(chain_a[:middle_chain_a], chain_b, similarity_func, gap_penalty)
        row_right = __nw_score_last_line(chain_a[middle_chain_a:][::-1], chain_b[::-1], similarity_func, gap_penalty)

        reversed_row_right = row_right[::-1]

        # Getting maximum
        row = [l + r for l, r in zip(row_left, reversed_row_right)]
        maxidx, maxval = max(enumerate(row), key=lambda a: a[1])

        middle_chain_b = maxidx

        # Recursive calls

        aligned_first_left, aligned_second_left = __hirschberg(chain_a[:middle_chain_a],
                                                               chain_b[:middle_chain_b], similarity_func, gap_penalty)
        algined_first_right, aligned_second_right = __hirschberg(chain_a[middle_chain_a:],
                                                                 chain_b[middle_chain_b:], similarity_func, gap_penalty)

        first_aligned = list(aligned_first_left) + list(algined_first_right)
        second_aligned = list(aligned_second_left) + list(aligned_second_right)

    return "".join(first_aligned), "".join(second_aligned)


def hirschberg(chain_a, chain_b, similarity_func, gap_penalty):
    return __hirschberg(chain_a, chain_b, similarity_func, gap_penalty)
