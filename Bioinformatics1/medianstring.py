'''
Median String Solution
'''

'''
This algorithm is designed to find a motif that is maximally conserved in all strands of a collection of DNA strands.
To accomplish this, the algorithm first generated all 4**k possible kmers. These kmers are then compared to an 
iteratively sliding window through each DNA strand to determine its degree of conservations. As each kmer is generated, 
it is only kept if it improves upon the previously best scoring kmer or is equally conserved. The resulting list 
of kmers represent motifs that may not be found exactly within any DNA strand, but are found with the least number 
of nucleotide changes possible.
'''



def hammingdistance(pattern, string):
    dist = 0
    for i in range(len(pattern)):
        if pattern[i] != string[i]:
            dist += 1
    return dist



def patternstringdistance(pattern, dna):
    k = len(pattern)
    distance = 0
    for seq in dna:
        l = len(seq)
        hamm_dist = float('inf')
        for i in range(l - k + 1):
            hamm_dist_window = hammingdistance(pattern, seq[i:i+k])
            if hamm_dist > hamm_dist_window:
                hamm_dist = hamm_dist_window
        distance += hamm_dist
    return distance



def numbertopattern(n, k):
    p = []
    seq1 = '0123ACGT'
    seq_dict = { int(seq1[i]):seq1[i + 4] for i in range(4) }
    for i in range(k):
        p.insert(0, seq_dict[n % 4])
        n //= 4
    return ''.join(p)



def allmedianstring(dna, k):
    dist = float('inf')
    for i in range(4**k):
        pattern = numbertopattern(i, k)
        curr_dist = patternstringdistance(pattern, dna)
        if dist > curr_dist:
            dist = curr_dist
            median = []
            median.append(pattern)
        elif dist == curr_dist:
            median.append(pattern)
    return median


'''
Example
'''

print(allmedianstring(['ACGCGCGCTT', 'ACGCGCGAGT', 'ACTGGCGCTT', 'CCGTGCGCAT'], 5))
# Output: ['CGCGC', 'GCGCG', 'GCGCT']
