'''
Most Probable kmer Solution
'''

'''
This algorithm determines the likely kmer within a given DNA sequence based on a collection of smaller sequences 
and a given length, k. A profile matrix is first constructed, which contains the relative frequency of each 
nucleotides at each position within the collection of smaller strands. This matrix is used to determine the
probability of every kmer within the larger DNA sequence. A kmer is only kept if it hsa a higher probability than
all previously tested kmers, so the final returned kmer represents the most probable kmer found within the DNA
sequence. It should be noted that if a nucleotide within a kmer of the DNA sequence does not occur at that
position, the kmers probability will be 0 even if it is highly likely to occur otherwise. Therefore, an alternate
method should be used to account for theoretically possible kmers.
'''



def profile(motifs):
    k = len(motifs[0])
    n = len(motifs)
    freq = 1 / n
    seq1 = 'ACGT0123'
    seq_dict = { seq1[i]:int(seq1[i+4]) for i in range(4) }
    matrix = [[0 for _ in range(k)] for _ in range(4)]
    for motif in motifs:
        for i in range(k):
            matrix[seq_dict[motif[i]]][i] += freq
    return matrix



def probability(pattern, profile):
    seq1 = 'ACGT0123'
    seq_dict = { seq1[i]:int(seq1[i+4]) for i in range(4) }
    prob = 1
    k = len(pattern)
    for i in range(k):
        prob *= profile[seq_dict[pattern[i]]][i]
    return prob



def profile_mprob_kmer(seq, k, motifs):
    prof = profile(motifs)
    l = len(seq)
    pmax = -float('inf')
    imax = -float('inf')
    for i in range(l - k + 1):
        prob = probability(seq[i:i + k], prof)
        if prob > pmax:
            pmax = prob
            imax = i
    return 'The most probable kmer is: ' + seq[imax:imax + k] + '\n' + 'Its probability is: ' + str(pmax)


'''
Example
'''

print(profile_mprob_kmer('TAACGTGCACTGATGCAG', 5, ['TTCAG', 'TAACT', 'ACCAG', 'CCCAG']))
# Output: The most probable kmer is: TAACG
# Its probability is: 0.005859375

# Note: 'TGCAG' is not returned because 'G' does not occur at the second nucleotide position.
