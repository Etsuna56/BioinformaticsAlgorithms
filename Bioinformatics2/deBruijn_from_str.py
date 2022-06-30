'''
deBruijn Graph from String Solution
'''

'''
A deBruijn graph is a representation of a overlapping segments among a collection of kmers where
nodes represent the overlapping prefixes and suffixes of kmers and edges joining two nodes represent
kmers formed by the overlapping prefix and suffix. This algorithm first forms all possible kmers of a given
nucleotide sequence, of which all kmer prefixes are added to a dictionary where keys are the graphical nodes 
that edges leave from and values are the nodes which edges end. These nodes and their edges are then
formatting with arrows to easily see the direction of each connections.
'''



def deBruijn_from_str(text, k): # creates dict of all nodes with their edges to other nodes
    l = len(text)
    adj = {}
    for i in range(l - k + 1):
        if text[i:i + k - 1] in adj:
            adj[text[i:i + k - 1]].append(text[i + 1:i + k]) # adds next kmer's prefix to previous kmer's suffix
        else:
            adj[text[i:i + k - 1]] = []
            adj[text[i:i + k - 1]].append(text[i + 1:i + k]) # adds next kmer's prefix to previous kmer's suffix
    connections = []
    for node1, node2 in adj.items():
        connections.append(node1 + ' -> ' + ' '.join(node2))
    deBruijn = '\n'.join(connections)
    return deBruijn



'''
Example
'''


print(deBruijn_from_str('ATGGATACG', 3))

# Output: AT -> TG TA
        # TG -> GG
        # GG -> GA
        # GA -> AT
        # TA -> AC
        # AC -> CG
