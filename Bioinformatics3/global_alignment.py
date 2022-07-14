'''
Global Alignment Solution
'''

'''
The following algorithm performs a global alignment between two strings by determining the longest path
in an alignment graph of the two strings. In this graph, edges are weighted to negatively score gaps and
mismatches in the alignment while positively valuing matches between the strings. As the maximum score of
each node is recorded in the scoring matrix, the edge used to achieve each maximum score is recorded in the
backtrack matrix for use in constructing the alignment of the strings.
'''


import numpy as np


def global_lcs_backtrack(v, w, sigma, mu): # v = vertical, w = horizontal
    n = len(v)
    m = len(w)

    s = np.matrix(np.zeros((n+1)*(m+1), dtype=int).reshape((n+1, m+1)))
    backtrack = np.matrix(np.zeros((n+1)*(m+1), dtype=int).reshape((n+1, m+1)))
    s[0] = [i*-sigma for i in range(m+1)]
    s[:, 0] = np.matrix([i*-sigma for i in range(n+1)]).reshape(n+1, 1)

    for i in range(1, n+1):
        for j in range(1, m+1):
            s[i, j] = max([s[i-1, j] - sigma, s[i, j-1] - sigma, s[i-1, j-1] + 1 if v[i-1] == w[j-1] else s[i-1, j-1] - mu]) # determines optimal path to current node from all preceeding nodes
            if s[i, j] == s[i-1, j] - sigma:
                backtrack[i, j] = 1 # 1 == down
            elif s[i, j] == s[i, j-1] - sigma:
                backtrack[i, j] = 2 # 2 == right
            elif s[i, j] == s[i-1, j-1] + 1 and v[i-1] == w[j-1]:
                backtrack[i, j] = 3 # 3 == diagonal match
            elif s[i, j] == s[i-1, j-1] - mu:
                backtrack[i, j] = 4 # 4 == diagonal mismatch

    return backtrack



def global_lcs_full(v, w, i, j, sigma, mu): # v = vertical, w = horizontal
    global lcs_v, lcs_w, score, backtrack

    backtrack = global_lcs_backtrack(v, w, sigma, mu)
    lcs_v = ''
    lcs_w = ''
    score = 0

    if i == 0 or j == 0:
        return

    def global_lcs_build(v, i, j):
        global lcs_v, lcs_w, score, backtrack

        while i > 0 and j > 0:
            if 1 == backtrack[i, j]: # backtracks up a single position
                i -= 1
                lcs_v += v[i]
                lcs_w += '-'
                score -= sigma
                continue
            elif 2 == backtrack[i, j]: # backtracks left a single position
                j -= 1
                lcs_w += w[j]
                lcs_v += '-'
                score -= sigma
                continue
            elif 3 == backtrack[i, j]: # backtracks diagonally, match added to lcs
                i -= 1
                j -= 1
                lcs_v += v[i]
                lcs_w += w[j]
                score += 1
                continue
            elif 4 == backtrack[i, j]:
                i -= 1
                j -= 1
                lcs_v += v[i]
                lcs_w += w[j]
                score -= mu
                continue

    if 1 == backtrack[i, j]: # determines first backtracking step to take
        score -= sigma
        lcs_v += v[i-1]
        lcs_w += '-'
        global_lcs_build(v, i - 1, j)
    elif 2 == backtrack[i, j]:
        score -= sigma
        lcs_w += w[j-1]
        lcs_v += '-'
        global_lcs_build(v, i, j - 1)
    elif 3 == backtrack[i, j]:
        score += 1
        lcs_v += v[i-1]
        lcs_w += w[j-1]
        global_lcs_build(v, i - 1, j - 1)
    elif 4 == backtrack[i, j]:
        score -= mu
        lcs_v += v[i-1]
        lcs_w += w[j-1]
        global_lcs_build(v, i - 1, j - 1)

    print('Alignment Score: ' + str(score))
    print(lcs_v[::-1])
    print(lcs_w[::-1])


'''
Example
'''


v = 'ATAGCGACGCCT'
w = 'ATACGATACA'

sigma = 1
mu = 2

global_lcs_full(v, w, len(v), len(w), sigma, mu)



# Output: Alignment Score: -1
# ATAGCGA--C-GCCT
# ATA-CGATACA----