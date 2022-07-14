'''
Longest Common Subsequence Solution
'''

'''
The following algorithm solves the longest common subsequence between two DNA strings by finding the longest
path through a grid representing the alignment of these two DNA strings. All edges in this grid have a zero
weight unless the nucleotides of both strings match, in which case the diagonal edges at this position yields
a score of one. A backtracking matrix is built by recording which edge yielded the longest path to each node
and using these edges to reverse through the strings and build the longest common subsequence.
'''


import numpy as np


def lcs_backtrack(v, w):
    n = len(v)
    m = len(w)

    s = np.matrix(np.zeros((n+1)*(m+1), dtype=int).reshape((n+1, m+1)))
    backtrack = np.matrix(np.zeros((n+1)*(m+1), dtype=int).reshape((n+1, m+1)))

    for i in range(1, n+1):
        for j in range(1, m+1):
            s[i, j] = max([s[i-1, j], s[i, j-1], s[i-1, j-1] + 1 if v[i-1] == w[j-1] else s[i-1, j-1]])
            if s[i, j] == s[i-1, j]:
                backtrack[i, j] = 1 # 1 == down
            elif s[i, j] == s[i, j-1]:
                backtrack[i, j] = 2 # 2 == right
            elif s[i, j] == s[i-1, j-1] + 1 and v[i-1] == w[j-1]:
                backtrack[i, j] = 3 # 3 == diagonal

    print('Scoring Matrix \n' + str(s))
    print('Backtracking Matrix: \n' + str(backtrack))
    return backtrack



def output_lcs(backtrack, v, i, j):

    lcs = ''
    if i == 0 or j == 0:
        return

    def lcs_build(backtrack, v, i, j):
        lcs = ''
        while i > 0 and j > 0:
            if 1 == backtrack[i, j]: # backtracks up a single position
                i -= 1
                continue
            elif 2 == backtrack[i, j]: # backtracks left a single position
                j -= 1
                continue
            else: # backtracks diagonally, match added to lcs
                i -= 1
                j -= 1
                lcs += v[i]
        return lcs

    if 1 == backtrack[i, j]: # determines first backtracking step to take
        lcs += lcs_build(backtrack, v, i - 1, j)
    elif 2 == backtrack[i, j]:
        lcs += lcs_build(backtrack, v, i, j - 1)
    else:
        lcs += lcs_build(backtrack, v, i - 1, j - 1)
        lcs += v[i-1]

    return 'Longest Common Subsequence: ' + lcs[::-1]


'''
Example
'''

v = 'TGTACG'
w = 'GCTAGT'

backtrack = lcs_backtrack(v, w)
print(output_lcs(backtrack, v, len(v), len(v)))


# Output: Scoring Matrix
# [[0 0 0 0 0 0 0]
#  [0 0 0 1 1 1 1]
#  [0 1 1 1 1 2 2]
#  [0 1 1 2 2 2 3]
#  [0 1 1 2 3 3 3]
#  [0 1 2 2 3 3 3]
#  [0 1 2 2 3 4 4]]

# Backtracking Matrix:
# [[0 0 0 0 0 0 0]
#  [0 1 1 3 2 2 2]
#  [0 3 2 1 1 3 2]
#  [0 1 1 3 2 1 3]
#  [0 1 1 1 3 2 1]
#  [0 1 3 1 1 1 1]
#  [0 1 1 1 1 3 2]]

# Longest Common Subsequence: GTAG
