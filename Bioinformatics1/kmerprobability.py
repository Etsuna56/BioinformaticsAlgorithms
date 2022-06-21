'''
Random Kmer Probability Solution
'''

'''
This algorithm calculates the probability of a motif of length, k, occurring
within a specified number of DNA strands, n, of a specific length, l.
'''



def kmerprobability(l, k, n):
    prob = float((0.25**k)*(l - k + 1)*(n))
    print('We expect ' + str(prob) + ' occurrences of the kmer motif within the DNA strands')


'''
Example
'''

kmerprobability(1000, 12, 500)
# Output: 'We expect 0.029474496841430664 occurrences of the kmer motif within the DNA strands'
