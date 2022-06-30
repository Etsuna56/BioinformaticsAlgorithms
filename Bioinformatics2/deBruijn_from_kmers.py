'''
deBruijn Graph from Kmers Solution
'''

'''
A deBruijn graph is a representation of a overlapping segments among a collection of kmers where
nodes represent the overlapping prefixes and suffixes of kmers and edges joining two nodes represent
kmers formed by the overlapping prefix and suffix. This algorithm forms a dictionary where keys are 
kmer prefixes and values are kmer suffixes. These keys and values are then formatted to easily
understand the directionality of all edges.
'''




def deBruijn_from_kmers(patterns):
    k = len(patterns[0])
    adj = {}
    for p in patterns:
        if p[:k - 1] in adj:
            adj[p[:k - 1]].append(p[1:])
        else:
            adj[p[:k - 1]] = []
            adj[p[:k - 1]].append(p[1:])
    connections = []
    for node1, node2 in adj.items():
        connections.append(node1 + ' -> ' + ' '.join(node2))
    deBruijn = '\n'.join(connections)
    return deBruijn



'''
Example
'''


print(deBruijn_from_kmers(['ATG', 'GCC', 'TGC', 'CCC', 'CCT']))

# Output: AT -> TG
        # GC -> CC
        # TG -> GC
        # CC -> CC CT
