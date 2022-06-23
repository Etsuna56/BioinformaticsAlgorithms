'''
Randomized Motif Search Solution
'''

'''
This set of functions is a Monte Carlo algorithm used to determine the most probable set of kmers within a collection
of DNA strands. A random set of kmers from each DNA strand is generated and a profile matrix is formed based on this
selection of kmers utilizing Laplace's rule to account for all theoretically possible kmers. This profile matrix
is then used to generate a new set of most probable kmers from each DNA strand, and this process iterates to form
increasingly more conserved sets of motifs. Due to its random formation of initial kmer sets, this algorithm must be
run a large amount of time to obtain better performance.
'''

import random



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



def biased_random(bias_list):
    number = random.uniform(0, sum(bias_list))
    curr = 0
    for i, bias in enumerate(bias_list):
        curr += bias
        if number <= curr:
            return i



def profile_rand_kmer(seq, k, profile):
    l = len(seq)
    p = []
    for i in range(l - k + 1):
        p.append(probability(seq[i:i+k], profile))
    i = biased_random(p)
    return seq[i:i + k]



def motifs_search(dna, k, profile):
    motifs = []
    for seq in dna:
        motifs.append(profile_mprob_kmer(seq, k, profile))
    return motifs



def one_rand_motifsearch(dna, k):
    t = len(dna)
    l = len(dna[0])
    m = [random.randint(0, l - k) for _ in range(t)]
    motifs = [dna[i][m[i]:m[i]+k] for i in range(t)]
    best_motifs = motifs
    best_score = score(best_motifs)
    while True:
        prof = pseudoprofile(motifs)
        motifs = motifs_search(dna, k, prof)
        curr_score = score(motifs)
        if curr_score < best_score:
            best_motifs = motifs
            best_score = curr_score
        else:
            return best_motifs, best_score



def rand_motifsearch(dna, k, _iter = 1000):
    t = len(dna)
    best_score = float('inf')
    random.seed()
    for _ in range(_iter):
        curr_bestmotifs, curr_bestscore = one_rand_motifsearch(dna, k)
        if curr_bestscore < best_score:
            best_motifs, best_score = (curr_bestmotifs, curr_bestscore)
    result = ''
    for motif in best_motifs:
        result += motif + ' '
    return result



'''
Example
'''

print(rand_motifsearch(['AACAGACGACAC', 'TGGAGACGAGGG', 'TTCCTACCCCAC', 'TTCTTTCGACAC'], 5))
# Output with Random Motif Search: 'GACAC GAGAC CCCAC GACAC' 
# Output with Greeedy Motif Search: 'GACAC GAGAC TACCC GACAC'

