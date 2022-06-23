'''
Distance Between String and Pattern Solution
'''

'''
Given a collection of DNA strands and a median string from the collection of kmers, this algorithm determines
the minimum number of nucleotide differences between the median string and all of the DNA strands. This is done
by sliding a window of equal length to the median string through each individual DNA strand and comparing the
window to the median string to determine the number of mismatch nucleotides. The minimum number of mismatches is
kept and added to a running count of total nucleotide differences, which is returned once the algorithm has
iterated through all DNA strands.
'''



def hammingdistance(pattern, string):
    dist = 0
    for i in range(len(pattern)):
        if pattern[i] != string[i]:
            dist += 1
    return dist



def patt_str_dist(pattern, dna):
    l = len(dna[0])
    k = len(pattern)
    dist = 0
    patterns = []
    for gene in dna:
        min_patt = ''
        hamdist = float('inf')
        for i in range(l - k + 1):
            patt = gene[i:i + k]
            if hammingdistance(patt, pattern) < hamdist:
                hamdist = hammingdistance(patt, pattern)
                min_patt = patt
        patterns.append(min_patt)
        dist += hamdist
    return 'Nucleotide mismatches: ' + str(dist) + '\n' + 'Patterns: ' + str(patterns)



'''
Example
'''

print(patt_str_dist('ACGGG', ['AACAGACGACAC', 'TGGAGACGAGGG', 'TTCCTACCCCAC', 'TTCTTTCGACAC']))
# Output: Nucleotide mismatches: 9
# Patterns: ['ACAGA', 'ACGAG', 'ACCCC', 'TCGAC']

