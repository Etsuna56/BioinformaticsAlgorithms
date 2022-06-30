'''
Overlap Graph Solution
'''

'''
Given a collection of kmers, this function returns an adjacency list in the form of a dictionary.
The keys of this adjacency list correspond to the nodes of an graph representing the overlapping
connetions between all kmers in the collections. This values of each key correspond to the
graphs edges, which are all kmers which the node shares an overlapping portion. This adjacency
list is formed by comparing the suffix of every kmer in the collection to every other kmer's prefix
and adds all matches to the kmer's list of edges.
'''


def overlapgraph(patterns):
    n = len(patterns)
    k = len(patterns[0])
    adj = {}
    for kmer in patterns:
        if kmer not in adj:
            adj[kmer] = []
    for i in range(n):
        for j in range(n):
            if patterns[i][1:] == patterns[j][:k-1]:
                adj[patterns[i]].append(patterns[j])
    answer = []
    for kmer in sorted(adj.keys()):
        if len(adj[kmer]) > 0:
            overlaps = ''
            for string in sorted(adj[kmer]):
                if string not in overlaps:
                    overlaps += ' ' + string
                    answer.append(kmer + ':' + overlaps)
    return answer



'''
Example
'''

print(overlapgraph(['ATG', 'GTT', 'TGT', 'TGA', 'AGG', 'GTT', 'ACC', 'ACT', 'CTT', 'TGC', 'GCT']))

# Output: ['ACT: CTT', 'ATG: TGA', 'ATG: TGA TGC', 'ATG: TGA TGC TGT', 'GCT: CTT', 'TGC: GCT', 'TGT: GTT']

