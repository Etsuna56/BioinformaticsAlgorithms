'''
Motif Enumeration Solution
'''

'''
This algorithm is designed to find all patterns of length k with at most d mismatches found within all DNA
strands of a given sample. It functions by iteratively sliding a window through the first strand in the sample and
generating all kmers with at most d mismatches to the window. A similar window is then iteratively slid through
each additional strand and all generated kmers are compared to the window in each additional strand.
If the generated kmer appears in all subsequent strings with at most d mismatches, it is added to a final
list of patterns and returned.
'''



def hammingdistance(pattern, string):
    dist = 0
    for i in range(len(pattern)):
        if pattern[i] != string[i]:
            dist += 1
    return dist


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



def motif_enumeration(dna_sample, k, d):
    patterns = []
    gene = dna_sample[0]
    for i in range(len(gene) - k + 1):
        neighbors_list = neighbors(gene[i:i+k], d)
        for neighbor in neighbors_list:
            count = 0
            for genes in dna_sample:
                for j in range(len(genes) - k + 1):
                    if hammingdistance(neighbor, genes[j:j + k]) <= d:
                        count += 1
                        break
            if count == len(dna_sample):
                patterns.append(neighbor)
    result = ''
    for pattern in patterns:
        result += str(pattern) + ' '
    return result


'''
Example
'''


print(motif_enumeration(['AGACGCCACTAT', 'AGGTCCCATATC', 'ATATGCCACTGA', 'ATCTTCCACTGT'], 5, 1))
# Output: 'TGCCA CCCAC CCCCT CCAAT CCATT CCACA'
