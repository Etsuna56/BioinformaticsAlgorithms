'''
Greedy Motif Search With Pseudocounts Solution
'''

'''
This greedy algorithm is designed to determine the most probable kmer in each string of a collection
of DNA strands by first building an initial profile matrix of nucleotide frequencies for the first kmer
in the first DNA strand, which is then used to determine the most probable kmer in each additional strand.
This profile utilizes Laplace's rule to account for theoretically possible strings in the formation of
most probable kmers found in each additional strand of DNA. As a new kmer is selected, the profile matrix 
is updated to reflect the addition of the kmers nucleotides. This final collection of kmers from all strands 
is then compared to the first set of first kmers in all strands to determine of a lower score has been achieved. 
If so, the new set of kmers in kept and the process repeats for all possible kmers in the first strand. 
'''



def pseudoprofile(motifs): # profile matrix with Laplace's rule
    k = len(motifs[0])
    n = len(motifs)
    freq = 1 / (n + 4)
    seq1 = 'ACGT0123'
    seq_dict = { seq1[i]:int(seq1[i+4]) for i in range(4) } # constructs dictionary to keep tally of nucleotides
    prof = [[1 for _ in range(k)] for __ in range(4)]  # creates lists for nucleotide probabilities (Laplace --> 1)
    for motif in motifs:
        for i in range(k):
            prof[seq_dict[motif[i]]][i] += freq
    return prof



def probability(pattern, profile):
    seq1 = 'ACGT0123'
    seq_dict = { seq1[i]:int(seq1[i+4]) for i in range(4) }
    prob = 1
    k = len(pattern)
    for i in range(k):
        prob *= profile[seq_dict[pattern[i]]][i]
    return prob



def profile_mprob_kmer(seq, k, prof):
    l = len(seq)
    pmax = -float('inf')
    imax = -float('inf')
    for i in range(l - k + 1):
        prob = probability(seq[i:i + k], prof)
        if prob > pmax:
            pmax = prob
            imax = i
    return seq[imax:imax + k]



def score(motifs):
    k = len(motifs[0])
    n = len(motifs)
    seq1 = 'ACGT0123'
    seq_dict = { seq1[i]:int(seq1[i+4]) for i in range(4) }
    p = [[0 for _ in range(4)] for __ in range(k)]
    for motif in motifs:
        for i in range(k):
            p[i][seq_dict[motif[i]]] += 1
    sm = 0
    for i in range(k):
        sm += max(p[i])
    return n * k - sm



def pseudogreedymotifsearch(dna, k):
    t = len(dna)
    dna1 = dna[0]
    l = len(dna1)
    best_motifs = [dna[i][0:k] for i in range(t)] # first kmers in each string of dna
    best_score = float('inf')
    for i in range(l - k + 1): # iterates across each nucleotide in dna strings
        motifs = []
        motifs.append(dna1[i:i+k]) # adds i-th kmer in 1st dna string to motifs
        for j in range(1, t): # iterates through each dna string
            prof = pseudoprofile(motifs) # creates profile matrix for i-th kmer of first strand
            motifs.append(profile_mprob_kmer(dna[j], k, prof)) # adds most prob kmers from dna strings to motifs
        curr_score = score(motifs)
        if curr_score < best_score: # compares score of most probable kmers
            best_motifs = motifs # replaces motifs if score is lowered
            best_score = curr_score # replaces score of new motifs
    result = ''
    for item in best_motifs:
        result += str(item) + ' '
    return result


'''
Example
'''

print(pseudogreedymotifsearch(['AACAGACGACAC', 'TGGAGACGAGGG', 'TTCCTACCCCAC', 'TTCTTTCGACAC'], 5))
# Output with Pseudocounts: 'GACAC GAGAC TACCC GACAC'
# Output without Pseudocounts: 'GACGA GACGA TTCCT TTCGA'

