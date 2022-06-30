'''
Read Pair String Reconstruction Solution
'''

'''
During genome sequencing, the length of reads if an important metric to optimize the assembly
process. Read pairs are one method that biologists employ to obtain more information from a
collection of DNA reads. Given a collection of paired reads, this algorithm sequences the genome
from which the read pairs were created. This is done by first creating a list of all paired read
prefixes and suffixes. These affixes represent the nodes of a deBruijn graph with the paired reads
themselves representing the graph's edges. The pairs of prefixes and suffixes forming the beginning
and end of a graphical edge are then assigned unique numbers identifying them as unique nodes for
path formation. A starting position is determined by finding the nodes with more outward edges than
inward edges. From this starting node, a depth first search is performed to find a Eulerian path
that crosses all edges. This path is used to reconstruct the original genome by concatenating a k + d
length prefix string to a suffix string and the suffix of the last paired kmer.
'''



import collections
from random import randint
from copy import deepcopy

def prefix(read_pair):
    pref = ''
    f = int((len(read_pair) - 1)/2)
    pref += (read_pair[:f - 1])
    pref += '|'
    pref += (read_pair[f + 1:len(read_pair) - 1])
    return pref

def suffix(read_pair):
    suff = ''
    f = int((len(read_pair) - 1)/2)
    suff += (read_pair[1:f])
    suff += '|'
    suff += (read_pair[len(read_pair) - f + 1:])
    return suff

def read_pair_assembly(reads, k, d):
    edges = [] # Format: [edge, [outward node, inward node]]
    for kmer in reads:
        edges.append(kmer)
        edges.append([prefix(kmer), suffix(kmer)])

    prefixes = [] # all prefix pairs
    suffixes = [] # all suffix pairs

    for i in range(0, len(edges), 2):
        prefixes.append(edges[i+1][0])
        suffixes.append(edges[i+1][1])

    pairs = {} # Format: {Prefix: [Suffix]}
    for i in range(len(prefixes)):
        pairs.setdefault(prefixes[i], [])
        pairs[prefixes[i]].append(suffixes[i])

    pairs = collections.OrderedDict(sorted(pairs.items()))
    numbers = {}
    rev_numbers = {}
    i = 0
    for pair in pairs: # assign pairs to unique IDs
        numbers[pair] = str(i)
        rev_numbers[str(i)] = pair
        i += 1

    for pair in pairs:
        for suff in pairs[pair]:
            numbers.setdefault(suff, [])
        if numbers[suff] == []:
            numbers[suff] = str(i)
            rev_numbers[str(i)] = suff
            i += 1

    adj_list = {}
    circuit_max = 0
    for pref in pairs:
        for suff in pairs[pref]:
            adj_list.setdefault(numbers[pref], [])
            adj_list[numbers[pref]].append(numbers[suff])
            circuit_max += 1

    red_adj_list = {}
    red_adj_list = deepcopy(adj_list)

    start = {}
    for key in red_adj_list:
        start.setdefault(key, 0)
        start[key] += len(red_adj_list[key])
    end = {}
    for key in red_adj_list:
        for value in red_adj_list[key]:
            end.setdefault(value, 0)
            end[value] += 1
    for key in end:
        try:
            if start[key] != end[key]:
                if start[key] > end[key]:
                    start_node = key
                if start[key] < end[key]:
                    end_node = key
        except KeyError:
            end_node = key
    for key in start:
        try:
            if end[key] != start[key]:
                if end[key] < start[key]:
                    start_node = key
                if end[key] > start[key]:
                    end_node = key
        except KeyError:
            start_node = key

    if end_node not in red_adj_list:
        red_adj_list[end_node] = []

    start = start_node
    curr_vert = start_node
    stack = [] # vertex  removed as edges are exhausted
    circuit = []
    while len(circuit) != circuit_max:
        if len(red_adj_list[curr_vert]) > 0:
            stack.append(curr_vert)
            pick = randint(0, len(red_adj_list[curr_vert])-1)
            temp = deepcopy(curr_vert)
            curr_vert = red_adj_list[temp][pick]
            red_adj_list[temp].remove(curr_vert)
        else:
            circuit.append(curr_vert)
            curr_vert = stack[len(stack) - 1]
            stack.pop()

    circuit = [start] + circuit[::-1]
    corr_order = []
    for vert in circuit:
        corr_order.append(rev_numbers[vert])

    a = []
    b = []
    for kmer in corr_order:
        a.append(kmer[:k-1])
        b.append(kmer[k:])

    prefixstring = ''
    for kmer in a:
        prefixstring += kmer[0]
    suffixstring = ''
    for kmer in b:
        suffixstring += kmer[0]

    genome = ''
    genome += prefixstring[0:k+d]
    genome += suffixstring
    genome += b[len(b) - 1][1:] # Adds suffix of final kmer
    return 'Sequenced genome: ' + genome


'''
Example
'''


print(read_pair_assembly(['AAC|ATC', 'ACG|TCG', 'AGC|AAC', 'CAT|ACT', 'CCT|CGC', 'CGC|CGA', 'CTA|GCA', 'GCA|GAC', 'GCC|ACG', 'TAA|CAT'], 3, 2))


# Output: 'Sequenced genome: AGCCTAACGCATCGACT'