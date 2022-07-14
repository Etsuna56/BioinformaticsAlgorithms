'''
Shared Kmers Solution
'''

'''
Given two genomic sequences, a kmer is shared if either itself or its reverse complement are found within both genomes.
The following algorithm forms a dictionary of all possible kmers in the first genome and stores their corresponding
coordinate positions. All kmers with reverse complements of the second genome are then checked to see if they occur
in the first, at which point matches are added to a set of coordinates to ensure than no duplicates are returned.
These shared kmers can then be graphically represented to visualize synteny blocks of the two genomes.
'''


from collections import defaultdict
import matplotlib.pyplot as plt


def defdict():
    return []


def reverse_complement(x):
    complement = ''
    for i in x:
        if i == 'A':
            complement += 'T'
        elif i == 'T':
            complement += 'A'
        elif i == 'C':
            complement += 'G'
        elif i == 'G':
            complement += 'C'
    return complement[::-1]


def shared_kmers(k, v, w):
    dct = defaultdict(defdict)
    for i in range(len(w)-k+1): # dict of all kmers in w with positions
        dct[w[i:i+k]].append(i)

    matches = set()
    for i in range(len(v)-k+1): # all kmers of v
        match = v[i:i+k]
        rev_match = reverse_complement(match)
        for compare in (match, rev_match):
            if compare in dct.keys(): # check if kmer or reverse comp in w
                for m in dct[compare]:
                    matches.add((i, m)) # add position of kmer in v and w as a tuple to matches

    print('Number of Matches: ' + str(len(matches)))
    print('Match Coordinates: ' + str(matches))

    return matches


'''
Example
'''



# shared_kmers(3, 'TCAGTTGGCCTACAT', 'CCTACATGAGGTCTG')


# Output: Number of Matches: 9
# Match Coordinates: {(12, 4), (8, 8), (1, 12), (11, 3), (8, 0), (0, 6), (12, 5), (9, 1), (10, 2)}


'''
Bacterial Genome Plotting
'''

# with open("e_coli_genome.txt", 'r') as v:
#     v = str(v.read().strip())
#     print('E. Coli Genome nucleotide count: ' + str(len(v)))
#     with open("salmonella_genome.txt", 'r') as w:
#         w = str(w.read().strip())
#         print('Salmonella Genome nucleotide count: ' + str(len(w)))
#
#
#     points = []
#     for point in shared_kmers(30, v, w):
#         points.append(point)
#
#     plt.scatter(*zip(*points))
#     plt.show()





