'''
Affine Gap Penalty Alignment Solution
'''

'''
Edits of large DNA strands are more likely to occur in chunks instead of one nucleotide at a time. Therefore, to
accurately penalize alignments of two strings, one must account different penalties for initiating gaps and
extending gaps. To implement this efficiently, separate backtracking and scoring matrices must be constructed
for each direction of edges. Transitioning between these matrices are scored differently than travelling within them
to penalize long gaps less severely than single nucleotide gaps.
'''


import numpy as np


'''
PAM250 Scoring Matrix
'''

pam = {'A': {'A': 2, 'C': -2, 'E': 0, 'D': 0, 'G': 1, 'F': -3, 'I': -1, 'H': -1, 'K': -1, 'M': -1, 'L': -2, 'N': 0, 'Q': 0, 'P': 1, 'S': 1, 'R': -2, 'T': 1, 'W': -6, 'V': 0, 'Y': -3}, 'C': {'A': -2, 'C': 12, 'E': -5, 'D': -5, 'G': -3, 'F': -4, 'I': -2, 'H': -3, 'K': -5, 'M': -5, 'L': -6, 'N': -4, 'Q': -5, 'P': -3, 'S': 0, 'R': -4, 'T': -2, 'W': -8, 'V': -2, 'Y': 0}, 'E': {'A': 0, 'C': -5, 'E': 4, 'D': 3, 'G': 0, 'F': -5, 'I': -2, 'H': 1, 'K': 0, 'M': -2, 'L': -3, 'N': 1, 'Q': 2, 'P': -1, 'S': 0, 'R': -1, 'T': 0, 'W': -7, 'V': -2, 'Y': -4}, 'D': {'A': 0, 'C': -5, 'E': 3, 'D': 4, 'G': 1, 'F': -6, 'I': -2, 'H': 1, 'K': 0, 'M': -3, 'L': -4, 'N': 2, 'Q': 2, 'P': -1, 'S': 0, 'R': -1, 'T': 0, 'W': -7, 'V': -2, 'Y': -4}, 'G': {'A': 1, 'C': -3, 'E': 0, 'D': 1, 'G': 5, 'F': -5, 'I': -3, 'H': -2, 'K': -2, 'M': -3, 'L': -4, 'N': 0, 'Q': -1, 'P': 0, 'S': 1, 'R': -3, 'T': 0, 'W': -7, 'V': -1, 'Y': -5}, 'F': {'A': -3, 'C': -4, 'E': -5, 'D': -6, 'G': -5, 'F': 9, 'I': 1, 'H': -2, 'K': -5, 'M': 0, 'L': 2, 'N': -3, 'Q': -5, 'P': -5, 'S': -3, 'R': -4, 'T': -3, 'W': 0, 'V': -1, 'Y': 7}, 'I': {'A': -1, 'C': -2, 'E': -2, 'D': -2, 'G': -3, 'F': 1, 'I': 5, 'H': -2, 'K': -2, 'M': 2, 'L': 2, 'N': -2, 'Q': -2, 'P': -2, 'S': -1, 'R': -2, 'T': 0, 'W': -5, 'V': 4, 'Y': -1}, 'H': {'A': -1, 'C': -3, 'E': 1, 'D': 1, 'G': -2, 'F': -2, 'I': -2, 'H': 6, 'K': 0, 'M': -2, 'L': -2, 'N': 2, 'Q': 3, 'P': 0, 'S': -1, 'R': 2, 'T': -1, 'W': -3, 'V': -2, 'Y': 0}, 'K': {'A': -1, 'C': -5, 'E': 0, 'D': 0, 'G': -2, 'F': -5, 'I': -2, 'H': 0, 'K': 5, 'M': 0, 'L': -3, 'N': 1, 'Q': 1, 'P': -1, 'S': 0, 'R': 3, 'T': 0, 'W': -3, 'V': -2, 'Y': -4}, 'M': {'A': -1, 'C': -5, 'E': -2, 'D': -3, 'G': -3, 'F': 0, 'I': 2, 'H': -2, 'K': 0, 'M': 6, 'L': 4, 'N': -2, 'Q': -1, 'P': -2, 'S': -2, 'R': 0, 'T': -1, 'W': -4, 'V': 2, 'Y': -2}, 'L': {'A': -2, 'C': -6, 'E': -3, 'D': -4, 'G': -4, 'F': 2, 'I': 2, 'H': -2, 'K': -3, 'M': 4, 'L': 6, 'N': -3, 'Q': -2, 'P': -3, 'S': -3, 'R': -3, 'T': -2, 'W': -2, 'V': 2, 'Y': -1}, 'N': {'A': 0, 'C': -4, 'E': 1, 'D': 2, 'G': 0, 'F': -3, 'I': -2, 'H': 2, 'K': 1, 'M': -2, 'L': -3, 'N': 2, 'Q': 1, 'P': 0, 'S': 1, 'R': 0, 'T': 0, 'W': -4, 'V': -2, 'Y': -2}, 'Q': {'A': 0, 'C': -5, 'E': 2, 'D': 2, 'G': -1, 'F': -5, 'I': -2, 'H': 3, 'K': 1, 'M': -1, 'L': -2, 'N': 1, 'Q': 4, 'P': 0, 'S': -1, 'R': 1, 'T': -1, 'W': -5, 'V': -2, 'Y': -4}, 'P': {'A': 1, 'C': -3, 'E': -1, 'D': -1, 'G': 0, 'F': -5, 'I': -2, 'H': 0, 'K': -1, 'M': -2, 'L': -3, 'N': 0, 'Q': 0, 'P': 6, 'S': 1, 'R': 0, 'T': 0, 'W': -6, 'V': -1, 'Y': -5}, 'S': {'A': 1, 'C': 0, 'E': 0, 'D': 0, 'G': 1, 'F': -3, 'I': -1, 'H': -1, 'K': 0, 'M': -2, 'L': -3, 'N': 1, 'Q': -1, 'P': 1, 'S': 2, 'R': 0, 'T': 1, 'W': -2, 'V': -1, 'Y': -3}, 'R': {'A': -2, 'C': -4, 'E': -1, 'D': -1, 'G': -3, 'F': -4, 'I': -2, 'H': 2, 'K': 3, 'M': 0, 'L': -3, 'N': 0, 'Q': 1, 'P': 0, 'S': 0, 'R': 6, 'T': -1, 'W': 2, 'V': -2, 'Y': -4}, 'T': {'A': 1, 'C': -2, 'E': 0, 'D': 0, 'G': 0, 'F': -3, 'I': 0, 'H': -1, 'K': 0, 'M': -1, 'L': -2, 'N': 0, 'Q': -1, 'P': 0, 'S': 1, 'R': -1, 'T': 3, 'W': -5, 'V': 0, 'Y': -3}, 'W': {'A': -6, 'C': -8, 'E': -7, 'D': -7, 'G': -7, 'F': 0, 'I': -5, 'H': -3, 'K': -3, 'M': -4, 'L': -2, 'N': -4, 'Q': -5, 'P': -6, 'S': -2, 'R': 2, 'T': -5, 'W': 17, 'V': -6, 'Y': 0}, 'V': {'A': 0, 'C': -2, 'E': -2, 'D': -2, 'G': -1, 'F': -1, 'I': 4, 'H': -2, 'K': -2, 'M': 2, 'L': 2, 'N': -2, 'Q': -2, 'P': -1, 'S': -1, 'R': -2, 'T': 0, 'W': -6, 'V': 4, 'Y': -2}, 'Y': {'A': -3, 'C': 0, 'E': -4, 'D': -4, 'G': -5, 'F': 7, 'I': -1, 'H': 0, 'K': -4, 'M': -2, 'L': -1, 'N': -2, 'Q': -4, 'P': -5, 'S': -3, 'R': -4, 'T': -3, 'W': 0, 'V': -2, 'Y': 10}}

'''
BLOSUM62 Scoring Matrix
'''

blosum = {'A': {'A': 4, 'C': 0, 'E': -1, 'D': -2, 'G': 0, 'F': -2, 'I': -1, 'H': -2, 'K': -1, 'M': -1, 'L': -1, 'N': -2, 'Q': -1, 'P': -1, 'S': 1, 'R': -1, 'T': 0, 'W': -3, 'V': 0, 'Y': -2}, 'C': {'A': 0, 'C': 9, 'E': -4, 'D': -3, 'G': -3, 'F': -2, 'I': -1, 'H': -3, 'K': -3, 'M': -1, 'L': -1, 'N': -3, 'Q': -3, 'P': -3, 'S': -1, 'R': -3, 'T': -1, 'W': -2, 'V': -1, 'Y': -2}, 'E': {'A': -1, 'C': -4, 'E': 5, 'D': 2, 'G': -2, 'F': -3, 'I': -3, 'H': 0, 'K': 1, 'M': -2, 'L': -3, 'N': 0, 'Q': 2, 'P': -1, 'S': 0, 'R': 0, 'T': -1, 'W': -3, 'V': -2, 'Y': -2}, 'D': {'A': -2, 'C': -3, 'E': 2, 'D': 6, 'G': -1, 'F': -3, 'I': -3, 'H': -1, 'K': -1, 'M': -3, 'L': -4, 'N': 1, 'Q': 0, 'P': -1, 'S': 0, 'R': -2, 'T': -1, 'W': -4, 'V': -3, 'Y': -3}, 'G': {'A': 0, 'C': -3, 'E': -2, 'D': -1, 'G': 6, 'F': -3, 'I': -4, 'H': -2, 'K': -2, 'M': -3, 'L': -4, 'N': 0, 'Q': -2, 'P': -2, 'S': 0, 'R': -2, 'T': -2, 'W': -2, 'V': -3, 'Y': -3}, 'F': {'A': -2, 'C': -2, 'E': -3, 'D': -3, 'G': -3, 'F': 6, 'I': 0, 'H': -1, 'K': -3, 'M': 0, 'L': 0, 'N': -3, 'Q': -3, 'P': -4, 'S': -2, 'R': -3, 'T': -2, 'W': 1, 'V': -1, 'Y': 3}, 'I': {'A': -1, 'C': -1, 'E': -3, 'D': -3, 'G': -4, 'F': 0, 'I': 4, 'H': -3, 'K': -3, 'M': 1, 'L': 2, 'N': -3, 'Q': -3, 'P': -3, 'S': -2, 'R': -3, 'T': -1, 'W': -3, 'V': 3, 'Y': -1}, 'H': {'A': -2, 'C': -3, 'E': 0, 'D': -1, 'G': -2, 'F': -1, 'I': -3, 'H': 8, 'K': -1, 'M': -2, 'L': -3, 'N': 1, 'Q': 0, 'P': -2, 'S': -1, 'R': 0, 'T': -2, 'W': -2, 'V': -3, 'Y': 2}, 'K': {'A': -1, 'C': -3, 'E': 1, 'D': -1, 'G': -2, 'F': -3, 'I': -3, 'H': -1, 'K': 5, 'M': -1, 'L': -2, 'N': 0, 'Q': 1, 'P': -1, 'S': 0, 'R': 2, 'T': -1, 'W': -3, 'V': -2, 'Y': -2}, 'M': {'A': -1, 'C': -1, 'E': -2, 'D': -3, 'G': -3, 'F': 0, 'I': 1, 'H': -2, 'K': -1, 'M': 5, 'L': 2, 'N': -2, 'Q': 0, 'P': -2, 'S': -1, 'R': -1, 'T': -1, 'W': -1, 'V': 1, 'Y': -1}, 'L': {'A': -1, 'C': -1, 'E': -3, 'D': -4, 'G': -4, 'F': 0, 'I': 2, 'H': -3, 'K': -2, 'M': 2, 'L': 4, 'N': -3, 'Q': -2, 'P': -3, 'S': -2, 'R': -2, 'T': -1, 'W': -2, 'V': 1, 'Y': -1}, 'N': {'A': -2, 'C': -3, 'E': 0, 'D': 1, 'G': 0, 'F': -3, 'I': -3, 'H': 1, 'K': 0, 'M': -2, 'L': -3, 'N': 6, 'Q': 0, 'P': -2, 'S': 1, 'R': 0, 'T': 0, 'W': -4, 'V': -3, 'Y': -2}, 'Q': {'A': -1, 'C': -3, 'E': 2, 'D': 0, 'G': -2, 'F': -3, 'I': -3, 'H': 0, 'K': 1, 'M': 0, 'L': -2, 'N': 0, 'Q': 5, 'P': -1, 'S': 0, 'R': 1, 'T': -1, 'W': -2, 'V': -2, 'Y': -1}, 'P': {'A': -1, 'C': -3, 'E': -1, 'D': -1, 'G': -2, 'F': -4, 'I': -3, 'H': -2, 'K': -1, 'M': -2, 'L': -3, 'N': -2, 'Q': -1, 'P': 7, 'S': -1, 'R': -2, 'T': -1, 'W': -4, 'V': -2, 'Y': -3}, 'S': {'A': 1, 'C': -1, 'E': 0, 'D': 0, 'G': 0, 'F': -2, 'I': -2, 'H': -1, 'K': 0, 'M': -1, 'L': -2, 'N': 1, 'Q': 0, 'P': -1, 'S': 4, 'R': -1, 'T': 1, 'W': -3, 'V': -2, 'Y': -2}, 'R': {'A': -1, 'C': -3, 'E': 0, 'D': -2, 'G': -2, 'F': -3, 'I': -3, 'H': 0, 'K': 2, 'M': -1, 'L': -2, 'N': 0, 'Q': 1, 'P': -2, 'S': -1, 'R': 5, 'T': -1, 'W': -3, 'V': -3, 'Y': -2}, 'T': {'A': 0, 'C': -1, 'E': -1, 'D': -1, 'G': -2, 'F': -2, 'I': -1, 'H': -2, 'K': -1, 'M': -1, 'L': -1, 'N': 0, 'Q': -1, 'P': -1, 'S': 1, 'R': -1, 'T': 5, 'W': -2, 'V': 0, 'Y': -2}, 'W': {'A': -3, 'C': -2, 'E': -3, 'D': -4, 'G': -2, 'F': 1, 'I': -3, 'H': -2, 'K': -3, 'M': -1, 'L': -2, 'N': -4, 'Q': -2, 'P': -4, 'S': -3, 'R': -3, 'T': -2, 'W': 11, 'V': -3, 'Y': 2}, 'V': {'A': 0, 'C': -1, 'E': -2, 'D': -3, 'G': -3, 'F': -1, 'I': 3, 'H': -3, 'K': -2, 'M': 1, 'L': 1, 'N': -3, 'Q': -2, 'P': -2, 'S': -2, 'R': -3, 'T': 0, 'W': -3, 'V': 4, 'Y': -1}, 'Y': {'A': -2, 'C': -2, 'E': -2, 'D': -3, 'G': -3, 'F': 3, 'I': -1, 'H': 2, 'K': -2, 'M': -1, 'L': -1, 'N': -2, 'Q': -1, 'P': -3, 'S': -2, 'R': -2, 'T': -2, 'W': 2, 'V': -1, 'Y': 7}}



def affine_gap_backtrack(v, w, match, mu, sigma, epsilon): # v = vertical, w = horizontal
    global backtrack, s_upper, s_middle, s_upper, b_lower, b_middle, b_upper

    n = len(v)
    m = len(w)

    s_upper = np.matrix(np.zeros((n+1)*(m+1), dtype=int).reshape((n+1, m+1)))
    s_middle = np.matrix(np.zeros((n+1)*(m+1), dtype=int).reshape((n+1, m+1)))
    s_lower = np.matrix(np.zeros((n+1)*(m+1), dtype=int).reshape((n+1, m+1)))

    b_upper = np.matrix(np.zeros((n+1)*(m+1), dtype=int).reshape((n+1, m+1)))
    b_middle = np.matrix(np.zeros((n + 1) * (m + 1), dtype=int).reshape((n + 1, m + 1)))
    b_lower = np.matrix(np.zeros((n + 1) * (m + 1), dtype=int).reshape((n + 1, m + 1)))


    s_lower[0] = [-10000 for i in range(m+1)]
    s_lower[:, 0] = np.matrix([-(sigma + epsilon*(i-1)) for i in range(n+1)]).reshape(n+1, 1)
    s_upper[:, 0] = np.matrix([-10000 for i in range(n+1)]).reshape(n+1, 1)
    s_upper[0] = [-(sigma + epsilon*(i-1)) for i in range(m+1)]
    s_middle[0] = [-(sigma + epsilon*(i-1)) for i in range(m+1)]
    s_middle[:, 0] = np.matrix([-(sigma + epsilon*(i-1)) for i in range(n+1)]).reshape(n+1, 1)
    s_lower[0, 0] = 0
    s_middle[0, 0] = 0
    s_upper[0, 0] = 0


    for i in range(1, n+1):
        for j in range(1, m+1):
            s_1 = s_lower[i-1, j] - epsilon # continue gap in w
            s_2 = s_middle[i-1, j] - sigma # start gap in w
            lower_score = max(s_1, s_2)
            s_lower[i, j] = lower_score
            if lower_score == s_2:
                b_lower[i, j] = 1 # initiate gap in w
            else:
                b_lower[i, j] = 0 # continue gap
            s_1 = s_upper[i, j-1] - epsilon # continue gap in v
            s_2 = s_middle[i, j-1] - sigma # start gap in v
            upper_score = max(s_1, s_2)
            s_upper[i, j] = upper_score
            if upper_score == s_2:
                b_upper[i, j] = 1 # initiate gap in v
            else:
                b_upper[i, j] = 0 # conmtinue gap
            s_3 = s_middle[i-1, j-1] + match if v[i-1] == w[j-1] else s_middle[i-1, j-1] - mu
            middle_score = max(lower_score, s_3, upper_score)
            s_middle[i, j] = middle_score
            if middle_score == lower_score:
                b_middle[i, j] = 0 # go to lower matrix
            elif middle_score == s_3:
                b_middle[i, j] = 1 # continue advancing in middle matrix
            elif middle_score == upper_score:
                b_middle[i, j] = 2 # go to upper matrix

    return s_middle[i, j]


def affine_gap_lcs(v, w, i, j, match, mu, sigma, epsilon): # v = vertical, w = horizontal
    global lcs_v, lcs_w, s_upper, s_middle, s_upper, b_lower, b_middle, b_upper

    score = affine_gap_backtrack(v, w, match, mu, sigma, epsilon)

    if i == 0 or j == 0:
        return

    def affine_lcs(v, i, j):
        global lcs_v, lcs_w, backtrack, s_upper, s_middle, s_upper, b_lower, b_middle, b_upper

        lcs_v = ''
        lcs_w = ''
        level = 1

        while i > 0 and j > 0:
            if level == 0:
                if b_lower[i, j] == 1: # initiate gap in w without extending then transition back to middle matrix
                    level = 1
                elif b_lower[i, j] == 0:
                    level = 0
                i -= 1
                lcs_v = v[i] + lcs_v
                lcs_w = '-' + lcs_w
                continue
            elif level == 2:
                if b_upper[i, j] == 1: # initiate gap in v without extending then back to middle matrix
                    level = 1
                elif b_upper[i, j] == 0:
                    level = 2
                j -= 1
                lcs_v = '-' + lcs_v
                lcs_w = w[j] + lcs_w
                continue
            else:
                if b_middle[i, j] == 1: # advance in middle matrix
                    i -= 1
                    j -= 1
                    lcs_v = v[i] + lcs_v
                    lcs_w = w[j] + lcs_w
                    continue
                else: # transition to upper or lower matrix
                    level = b_middle[i, j]

    affine_lcs(v, i, j)

    print('Score with Affine Gap Penalty: ' + str(score))
    print(lcs_v)
    print(lcs_w)
    

'''
Example
'''


v = 'PRTEINS'
w = 'PRTWPSEIN'

match = 6
mu = 2
sigma = 11
epsilon = 1

affine_gap_lcs(v, w, len(v), len(w), match, mu, sigma, epsilon)


# Score with Affine Gap Penalty: 12
# PRT---EINS
# PRTWPSEIN-