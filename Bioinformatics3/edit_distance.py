'''
Editing Distance Solution
'''

'''
The editing distance is equivalent to the minimum number of editing operation required to transform one string
into another. Therefore, in an alignment graph it is the total number of mismatches and indels, or gaps, in the
global alignment between the two strings. This is easily found by forming the scoring matrix where all mismatches
and indels are weighted with a score of 1 and all matches are zero weight edges.
'''


import numpy as np


'''
PAM250 Scoring Matrix
'''

pam = {'A': {'A': 2, 'C': -2, 'E': 0, 'D': 0, 'G': 1, 'F': -3, 'I': -1, 'H': -1, 'K': -1, 'M': -1, 'L': -2, 'N': 0, 'Q': 0, 'P': 1, 'S': 1, 'R': -2, 'T': 1, 'W': -6, 'V': 0, 'Y': -3}, 'C': {'A': -2, 'C': 12, 'E': -5, 'D': -5, 'G': -3, 'F': -4, 'I': -2, 'H': -3, 'K': -5, 'M': -5, 'L': -6, 'N': -4, 'Q': -5, 'P': -3, 'S': 0, 'R': -4, 'T': -2, 'W': -8, 'V': -2, 'Y': 0}, 'E': {'A': 0, 'C': -5, 'E': 4, 'D': 3, 'G': 0, 'F': -5, 'I': -2, 'H': 1, 'K': 0, 'M': -2, 'L': -3, 'N': 1, 'Q': 2, 'P': -1, 'S': 0, 'R': -1, 'T': 0, 'W': -7, 'V': -2, 'Y': -4}, 'D': {'A': 0, 'C': -5, 'E': 3, 'D': 4, 'G': 1, 'F': -6, 'I': -2, 'H': 1, 'K': 0, 'M': -3, 'L': -4, 'N': 2, 'Q': 2, 'P': -1, 'S': 0, 'R': -1, 'T': 0, 'W': -7, 'V': -2, 'Y': -4}, 'G': {'A': 1, 'C': -3, 'E': 0, 'D': 1, 'G': 5, 'F': -5, 'I': -3, 'H': -2, 'K': -2, 'M': -3, 'L': -4, 'N': 0, 'Q': -1, 'P': 0, 'S': 1, 'R': -3, 'T': 0, 'W': -7, 'V': -1, 'Y': -5}, 'F': {'A': -3, 'C': -4, 'E': -5, 'D': -6, 'G': -5, 'F': 9, 'I': 1, 'H': -2, 'K': -5, 'M': 0, 'L': 2, 'N': -3, 'Q': -5, 'P': -5, 'S': -3, 'R': -4, 'T': -3, 'W': 0, 'V': -1, 'Y': 7}, 'I': {'A': -1, 'C': -2, 'E': -2, 'D': -2, 'G': -3, 'F': 1, 'I': 5, 'H': -2, 'K': -2, 'M': 2, 'L': 2, 'N': -2, 'Q': -2, 'P': -2, 'S': -1, 'R': -2, 'T': 0, 'W': -5, 'V': 4, 'Y': -1}, 'H': {'A': -1, 'C': -3, 'E': 1, 'D': 1, 'G': -2, 'F': -2, 'I': -2, 'H': 6, 'K': 0, 'M': -2, 'L': -2, 'N': 2, 'Q': 3, 'P': 0, 'S': -1, 'R': 2, 'T': -1, 'W': -3, 'V': -2, 'Y': 0}, 'K': {'A': -1, 'C': -5, 'E': 0, 'D': 0, 'G': -2, 'F': -5, 'I': -2, 'H': 0, 'K': 5, 'M': 0, 'L': -3, 'N': 1, 'Q': 1, 'P': -1, 'S': 0, 'R': 3, 'T': 0, 'W': -3, 'V': -2, 'Y': -4}, 'M': {'A': -1, 'C': -5, 'E': -2, 'D': -3, 'G': -3, 'F': 0, 'I': 2, 'H': -2, 'K': 0, 'M': 6, 'L': 4, 'N': -2, 'Q': -1, 'P': -2, 'S': -2, 'R': 0, 'T': -1, 'W': -4, 'V': 2, 'Y': -2}, 'L': {'A': -2, 'C': -6, 'E': -3, 'D': -4, 'G': -4, 'F': 2, 'I': 2, 'H': -2, 'K': -3, 'M': 4, 'L': 6, 'N': -3, 'Q': -2, 'P': -3, 'S': -3, 'R': -3, 'T': -2, 'W': -2, 'V': 2, 'Y': -1}, 'N': {'A': 0, 'C': -4, 'E': 1, 'D': 2, 'G': 0, 'F': -3, 'I': -2, 'H': 2, 'K': 1, 'M': -2, 'L': -3, 'N': 2, 'Q': 1, 'P': 0, 'S': 1, 'R': 0, 'T': 0, 'W': -4, 'V': -2, 'Y': -2}, 'Q': {'A': 0, 'C': -5, 'E': 2, 'D': 2, 'G': -1, 'F': -5, 'I': -2, 'H': 3, 'K': 1, 'M': -1, 'L': -2, 'N': 1, 'Q': 4, 'P': 0, 'S': -1, 'R': 1, 'T': -1, 'W': -5, 'V': -2, 'Y': -4}, 'P': {'A': 1, 'C': -3, 'E': -1, 'D': -1, 'G': 0, 'F': -5, 'I': -2, 'H': 0, 'K': -1, 'M': -2, 'L': -3, 'N': 0, 'Q': 0, 'P': 6, 'S': 1, 'R': 0, 'T': 0, 'W': -6, 'V': -1, 'Y': -5}, 'S': {'A': 1, 'C': 0, 'E': 0, 'D': 0, 'G': 1, 'F': -3, 'I': -1, 'H': -1, 'K': 0, 'M': -2, 'L': -3, 'N': 1, 'Q': -1, 'P': 1, 'S': 2, 'R': 0, 'T': 1, 'W': -2, 'V': -1, 'Y': -3}, 'R': {'A': -2, 'C': -4, 'E': -1, 'D': -1, 'G': -3, 'F': -4, 'I': -2, 'H': 2, 'K': 3, 'M': 0, 'L': -3, 'N': 0, 'Q': 1, 'P': 0, 'S': 0, 'R': 6, 'T': -1, 'W': 2, 'V': -2, 'Y': -4}, 'T': {'A': 1, 'C': -2, 'E': 0, 'D': 0, 'G': 0, 'F': -3, 'I': 0, 'H': -1, 'K': 0, 'M': -1, 'L': -2, 'N': 0, 'Q': -1, 'P': 0, 'S': 1, 'R': -1, 'T': 3, 'W': -5, 'V': 0, 'Y': -3}, 'W': {'A': -6, 'C': -8, 'E': -7, 'D': -7, 'G': -7, 'F': 0, 'I': -5, 'H': -3, 'K': -3, 'M': -4, 'L': -2, 'N': -4, 'Q': -5, 'P': -6, 'S': -2, 'R': 2, 'T': -5, 'W': 17, 'V': -6, 'Y': 0}, 'V': {'A': 0, 'C': -2, 'E': -2, 'D': -2, 'G': -1, 'F': -1, 'I': 4, 'H': -2, 'K': -2, 'M': 2, 'L': 2, 'N': -2, 'Q': -2, 'P': -1, 'S': -1, 'R': -2, 'T': 0, 'W': -6, 'V': 4, 'Y': -2}, 'Y': {'A': -3, 'C': 0, 'E': -4, 'D': -4, 'G': -5, 'F': 7, 'I': -1, 'H': 0, 'K': -4, 'M': -2, 'L': -1, 'N': -2, 'Q': -4, 'P': -5, 'S': -3, 'R': -4, 'T': -3, 'W': 0, 'V': -2, 'Y': 10}}



def edit_dist_backtrack(v, w, sigma, mu): # v = vertical, w = horizontal
    n = len(v)
    m = len(w)

    s = np.matrix(np.zeros((n+1)*(m+1), dtype=int).reshape((n+1, m+1)))
    backtrack = np.matrix(np.zeros((n+1)*(m+1), dtype=int).reshape((n+1, m+1)))
    s[0] = [i*-sigma for i in range(m+1)]
    s[:, 0] = np.matrix([i*-sigma for i in range(n+1)]).reshape(n+1, 1)

    for i in range(1, n+1):
        for j in range(1, m+1):
            s[i, j] = max([s[i-1, j] - sigma, s[i, j-1] - sigma, s[i-1, j-1] if v[i-1] == w[j-1] else s[i-1, j-1] - mu]) # determines optimal path to current node from all preceeding nodes
            if s[i, j] == s[i-1, j] - sigma:
                backtrack[i, j] = 1 # 1 == down
            elif s[i, j] == s[i, j-1] - sigma:
                backtrack[i, j] = 2 # 2 == right
            elif s[i, j] == s[i-1, j-1] and v[i-1] == w[j-1]:
                backtrack[i, j] = 3 # 3 == diagonal match
            elif s[i, j] == s[i-1, j-1] - mu:
                backtrack[i, j] = 4 # 4 == diagonal mismatch
    return backtrack



def editing_distance(v, w, i, j, sigma, mu): # v = vertical, w = horizontal
    global lcs_v, lcs_w, score, backtrack

    backtrack = edit_dist_backtrack(v, w, sigma, mu)
    lcs_v = ''
    lcs_w = ''
    score = 0

    if i == 0 or j == 0:
        return

    def global_lcs(v, i, j):
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
        global_lcs(v, i - 1, j)
    elif 2 == backtrack[i, j]:
        score -= sigma
        lcs_w += w[j-1]
        lcs_v += '-'
        global_lcs(v, i, j - 1)
    elif 3 == backtrack[i, j]:
        lcs_v += v[i-1]
        lcs_w += w[j-1]
        global_lcs(v, i - 1, j - 1)
    elif 4 == backtrack[i, j]:
        score -= mu
        lcs_v += v[i-1]
        lcs_w += w[j-1]
        global_lcs(v, i - 1, j - 1)

    print('Editing Distance: ' + str(abs(score)))
    print(lcs_v[::-1])
    print(lcs_w[::-1])

'''
Example
'''


v = 'GATGGACTC'
w = 'GAGAATACC'

sigma = 1
mu = 1

editing_distance(v, w, len(v), len(w), 1, 1)


# Output: Editing Distance: 5
# GATGGACT-C-
# GA-GAA-TACC