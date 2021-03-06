'''
Neighbor Generation Solution
'''

'''
This algorithm is designed to generate all possible strings with at most d mismatches from a given pattern.
It functions by first checking that a maximum number of mismatches is given and that the given pattern is longer
than a single nucleotide. A set of 'neighbors' is then generated by recursively taking a single nucleotide
from the original pattern and creating a list of strings with alternative nucleotides to the pattern at that
position.
'''



def neighbors(pattern, d):
    nucleotides = ['A', 'G', 'C', 'T']
    neighborhood = []
    if d == 0:
        return [pattern]
    if len(pattern) == 1:
        return nucleotides
    suff_neighbors = neighbors(pattern[1:], d-1)
    for suff_neighbor in suff_neighbors:
        for i in nucleotides:
            if i != pattern[0]:
                neighborhood += [i + suff_neighbor]
    if d < len(pattern):
        suff_neighbors = neighbors(pattern[1:], d)
        for suff_neighbor in suff_neighbors:
            neighborhood += [pattern[0] + suff_neighbor]
    if pattern not in neighborhood:
        neighborhood.append(pattern)
    return neighborhood



'''
Example
'''


print(neighbors('GCA', 1))
# ['ACA', 'CCA', 'TCA', 'GAA', 'GGA', 'GTA', 'GCA', 'GCG', 'GCC', 'GCT']
