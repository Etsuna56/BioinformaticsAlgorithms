'''
Entropy Score Calculator Solution
'''

'''
This algorithm calculates the uncertainty of the probability distribution of a collection 
of motifs by finding the relative frequency, x, of each nucleotide at each position of the 
motifs and plugging those frequencies into the following entropy formula:

                        entropy = sum(-x * log2(x))
'''


import math



def entropy(prob):
    return sum([-x*math.log2(x) for x in prob if x > 0])



def entropy_motif(motifs):
    entropy_score = 0
    for i in range(len(motifs[0])):
        freq = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
        for motif in motifs:
            freq[motif[i]] += 1
        entropy_score += entropy([x/len(motifs) for x in freq.values()]) # Plugs relative freq of each nucleotide into ent
    print('The entropy score is: ' + str(round(entropy_score, 6)))


'''
Example
'''

entropy_motif(['ACGCGC', 'CCCCCC', 'GTCAAA', 'ACGTTC'])
# Output: 'The entropy score is: 7.622556'
