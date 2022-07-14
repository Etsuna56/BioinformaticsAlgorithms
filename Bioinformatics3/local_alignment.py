'''
Local Alignment Solution
'''

'''
The following algorithm performs a local alignment between two strings by determining the longest path
in an alignment graph of the two strings. In this graph, edges are weighted to negatively score gaps and
mismatches in the alignment while positively valuing matches between the strings. However, every node has an
additional edge connecting it directly to the source of the graph to ensure segments of little similarity between
the two strings do not impact the score of highly conserved subsequences. As the maximum score of each node is 
recorded in the scoring matrix, the edge used to achieve each maximum score is recorded in the backtrack matrix 
for use in constructing the alignment of the strings. Backtracking is initiated from the node which achieved the
highest overall score in the entire graph, which represents the end of locally aligned subsequence of the two strings.
'''


import numpy as np


'''
PAM250 Scoring Matrix
'''

pam = {'A': {'A': 2, 'C': -2, 'E': 0, 'D': 0, 'G': 1, 'F': -3, 'I': -1, 'H': -1, 'K': -1, 'M': -1, 'L': -2, 'N': 0, 'Q': 0, 'P': 1, 'S': 1, 'R': -2, 'T': 1, 'W': -6, 'V': 0, 'Y': -3}, 'C': {'A': -2, 'C': 12, 'E': -5, 'D': -5, 'G': -3, 'F': -4, 'I': -2, 'H': -3, 'K': -5, 'M': -5, 'L': -6, 'N': -4, 'Q': -5, 'P': -3, 'S': 0, 'R': -4, 'T': -2, 'W': -8, 'V': -2, 'Y': 0}, 'E': {'A': 0, 'C': -5, 'E': 4, 'D': 3, 'G': 0, 'F': -5, 'I': -2, 'H': 1, 'K': 0, 'M': -2, 'L': -3, 'N': 1, 'Q': 2, 'P': -1, 'S': 0, 'R': -1, 'T': 0, 'W': -7, 'V': -2, 'Y': -4}, 'D': {'A': 0, 'C': -5, 'E': 3, 'D': 4, 'G': 1, 'F': -6, 'I': -2, 'H': 1, 'K': 0, 'M': -3, 'L': -4, 'N': 2, 'Q': 2, 'P': -1, 'S': 0, 'R': -1, 'T': 0, 'W': -7, 'V': -2, 'Y': -4}, 'G': {'A': 1, 'C': -3, 'E': 0, 'D': 1, 'G': 5, 'F': -5, 'I': -3, 'H': -2, 'K': -2, 'M': -3, 'L': -4, 'N': 0, 'Q': -1, 'P': 0, 'S': 1, 'R': -3, 'T': 0, 'W': -7, 'V': -1, 'Y': -5}, 'F': {'A': -3, 'C': -4, 'E': -5, 'D': -6, 'G': -5, 'F': 9, 'I': 1, 'H': -2, 'K': -5, 'M': 0, 'L': 2, 'N': -3, 'Q': -5, 'P': -5, 'S': -3, 'R': -4, 'T': -3, 'W': 0, 'V': -1, 'Y': 7}, 'I': {'A': -1, 'C': -2, 'E': -2, 'D': -2, 'G': -3, 'F': 1, 'I': 5, 'H': -2, 'K': -2, 'M': 2, 'L': 2, 'N': -2, 'Q': -2, 'P': -2, 'S': -1, 'R': -2, 'T': 0, 'W': -5, 'V': 4, 'Y': -1}, 'H': {'A': -1, 'C': -3, 'E': 1, 'D': 1, 'G': -2, 'F': -2, 'I': -2, 'H': 6, 'K': 0, 'M': -2, 'L': -2, 'N': 2, 'Q': 3, 'P': 0, 'S': -1, 'R': 2, 'T': -1, 'W': -3, 'V': -2, 'Y': 0}, 'K': {'A': -1, 'C': -5, 'E': 0, 'D': 0, 'G': -2, 'F': -5, 'I': -2, 'H': 0, 'K': 5, 'M': 0, 'L': -3, 'N': 1, 'Q': 1, 'P': -1, 'S': 0, 'R': 3, 'T': 0, 'W': -3, 'V': -2, 'Y': -4}, 'M': {'A': -1, 'C': -5, 'E': -2, 'D': -3, 'G': -3, 'F': 0, 'I': 2, 'H': -2, 'K': 0, 'M': 6, 'L': 4, 'N': -2, 'Q': -1, 'P': -2, 'S': -2, 'R': 0, 'T': -1, 'W': -4, 'V': 2, 'Y': -2}, 'L': {'A': -2, 'C': -6, 'E': -3, 'D': -4, 'G': -4, 'F': 2, 'I': 2, 'H': -2, 'K': -3, 'M': 4, 'L': 6, 'N': -3, 'Q': -2, 'P': -3, 'S': -3, 'R': -3, 'T': -2, 'W': -2, 'V': 2, 'Y': -1}, 'N': {'A': 0, 'C': -4, 'E': 1, 'D': 2, 'G': 0, 'F': -3, 'I': -2, 'H': 2, 'K': 1, 'M': -2, 'L': -3, 'N': 2, 'Q': 1, 'P': 0, 'S': 1, 'R': 0, 'T': 0, 'W': -4, 'V': -2, 'Y': -2}, 'Q': {'A': 0, 'C': -5, 'E': 2, 'D': 2, 'G': -1, 'F': -5, 'I': -2, 'H': 3, 'K': 1, 'M': -1, 'L': -2, 'N': 1, 'Q': 4, 'P': 0, 'S': -1, 'R': 1, 'T': -1, 'W': -5, 'V': -2, 'Y': -4}, 'P': {'A': 1, 'C': -3, 'E': -1, 'D': -1, 'G': 0, 'F': -5, 'I': -2, 'H': 0, 'K': -1, 'M': -2, 'L': -3, 'N': 0, 'Q': 0, 'P': 6, 'S': 1, 'R': 0, 'T': 0, 'W': -6, 'V': -1, 'Y': -5}, 'S': {'A': 1, 'C': 0, 'E': 0, 'D': 0, 'G': 1, 'F': -3, 'I': -1, 'H': -1, 'K': 0, 'M': -2, 'L': -3, 'N': 1, 'Q': -1, 'P': 1, 'S': 2, 'R': 0, 'T': 1, 'W': -2, 'V': -1, 'Y': -3}, 'R': {'A': -2, 'C': -4, 'E': -1, 'D': -1, 'G': -3, 'F': -4, 'I': -2, 'H': 2, 'K': 3, 'M': 0, 'L': -3, 'N': 0, 'Q': 1, 'P': 0, 'S': 0, 'R': 6, 'T': -1, 'W': 2, 'V': -2, 'Y': -4}, 'T': {'A': 1, 'C': -2, 'E': 0, 'D': 0, 'G': 0, 'F': -3, 'I': 0, 'H': -1, 'K': 0, 'M': -1, 'L': -2, 'N': 0, 'Q': -1, 'P': 0, 'S': 1, 'R': -1, 'T': 3, 'W': -5, 'V': 0, 'Y': -3}, 'W': {'A': -6, 'C': -8, 'E': -7, 'D': -7, 'G': -7, 'F': 0, 'I': -5, 'H': -3, 'K': -3, 'M': -4, 'L': -2, 'N': -4, 'Q': -5, 'P': -6, 'S': -2, 'R': 2, 'T': -5, 'W': 17, 'V': -6, 'Y': 0}, 'V': {'A': 0, 'C': -2, 'E': -2, 'D': -2, 'G': -1, 'F': -1, 'I': 4, 'H': -2, 'K': -2, 'M': 2, 'L': 2, 'N': -2, 'Q': -2, 'P': -1, 'S': -1, 'R': -2, 'T': 0, 'W': -6, 'V': 4, 'Y': -2}, 'Y': {'A': -3, 'C': 0, 'E': -4, 'D': -4, 'G': -5, 'F': 7, 'I': -1, 'H': 0, 'K': -4, 'M': -2, 'L': -1, 'N': -2, 'Q': -4, 'P': -5, 'S': -3, 'R': -4, 'T': -3, 'W': 0, 'V': -2, 'Y': 10}}



def local_lcs_backtrack(v, w, sigma, matrix):
    global s, n, m

    n = len(v)
    m = len(w)

    s = np.matrix(np.zeros((n+1)*(m+1), dtype=int).reshape((n+1, m+1)))
    backtrack = np.matrix(np.zeros((n+1)*(m+1), dtype=int).reshape((n+1, m+1)))

    for i in range(1, n+1):
        for j in range(1, m+1):
            s[i, j] = max([s[i - 1, j] - sigma, s[i, j - 1] - sigma, s[i - 1, j - 1] + matrix[v[i - 1]][w[j - 1]], 0])  # determines optimal path to current node from all preceding nodes
            if s[i, j] == s[i - 1, j] - sigma:
                backtrack[i, j] = 1  # 1 == down
            elif s[i, j] == s[i, j - 1] - sigma:
                backtrack[i, j] = 2  # 2 == right
            elif s[i, j] == s[i - 1, j - 1] + matrix[v[i - 1]][w[j - 1]]:
                backtrack[i, j] = 3  # 3 == diagonal match
            elif s[i, j] == 0:
                backtrack[i, j] = 4  # 4 == back to start node

    return backtrack



def local_lcs_full(v, w, i, j, sigma, matrix):
    global best_substrings, best_score, backtrack, substring_v, substring_w, score, s, n, m

    backtrack = local_lcs_backtrack(v, w, sigma, matrix)
    max_coords = np.where(s == np.amax(s))
    max_i, max_j = list(zip(max_coords[0], max_coords[1]))[0]
    best_substrings = []
    substring_v = ''
    substring_w = ''
    best_score = -float('inf')

    if i == 0 or j == 0:
        return


    def local_lcs(v, i, j, sigma, matrix):
        global best_substrings, best_score, backtrack, substring_v, substring_w, score

        substring_v = v[i]
        susbtring_w = w[j]

        while i > 0 and j > 0:
            if 1 == backtrack[i, j]: # backtracks up a single position
                i -= 1
                substring_v += v[i]
                substring_w += '-'
                score -= sigma
                continue
            elif 2 == backtrack[i, j]: # backtracks left a single position
                j -= 1
                substring_v += '-'
                substring_w += w[j]
                score -= sigma
                continue
            elif 3 == backtrack[i, j]: # backtracks diagonally, match added to lcs
                i -= 1
                j -= 1
                score += matrix[v[i]][w[j]]
                substring_v += v[i]
                substring_w += w[j]
                continue
            elif 4 == backtrack[i, j]:
                i = 0
                j = 0
        if score > best_score:
            best_score = score
            best_substrings = [substring_v[::-1] + ':' + substring_w[::-1]]



    if 1 == backtrack[max_i, max_j]: # determines first backtracking step to take
        score = -sigma
        substring_v += v[max_i-1]
        substring_w += '-'
        local_lcs(v, max_i - 1, max_j, sigma, matrix)
    elif 2 == backtrack[max_i, max_j]:
        score = -sigma
        substring_v += '-'
        substring_w += w[max_j-1]
        local_lcs(v, max_i, max_j - 1, sigma, matrix)
    elif 3 == backtrack[max_i, max_j]:
        score = matrix[v[max_i-1]][w[max_j-1]]
        substring_v += v[max_i-1]
        substring_w += w[max_j-1]
        local_lcs(v, max_i - 1, max_j - 1, sigma, matrix)

    print('Local Alignment Score: ' + str(best_score))
    print(best_substrings[0].split(':')[0])
    print(best_substrings[0].split(':')[1])


'''
Example
'''


v = 'MEANLY'
w = 'PENALTY'

sigma = 1


local_lcs_full(v, w, len(v), len(w), sigma, pam)


# Output: Local Alignment Score: 19
# E-ANL-Y
# ENA-LTY