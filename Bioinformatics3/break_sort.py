'''
2-Break Sorting Solution
'''

'''
The 2-Break Sorting Algorithm below reconstructs a collection of 2-break operation used to transform genome q into
genome p. To guarantee that each operation is optimal. The original collection of alternating cycles is formed. From
this collection of cycles, a non-trivial cycle, which represents a cycle with more than 2 edges, is chosen and two
edges from genome q are broken to form new edges, one of which is guaranteed to be a part of genome p, which will
form a new trivial cycle consisting of two equivalent edges. Once all cycles are trivial, genome q and p will be
equivalent and the sorting will be complete.
*Note that chromosomes within the genomes are assumed to by cyclic so rearrangements and genome outputs may seem
to differ from the correct order. However, this is simply due to an alternate circular ordering.
'''


def chrom_to_cycle(chromosome): # chromosomal permutation to graphical nodes
    l = len(chromosome)
    nodes = [0]*2*l
    for i in range(l):
        try:
            block = chromosome[i]
            num = int(block[1:])
            if block[0] == '+':
                nodes[2*i] = num*2 - 1
                nodes[2*i+1] = num*2
            elif block[0] == '-':
                nodes[2*i] = 2*num
                nodes[2*i+1] = 2*num - 1
        except ValueError:
            break
    # print(nodes)
    return nodes


def cycle_to_chrom(cycle): # graphical cycle to chromosomal permutation
    l = len(cycle)//2
    chromosome = [0]*l
    for i in range(l):
        if cycle[2*i] < cycle[2*i+1]:
            chromosome[i] = '+' + str(cycle[2*i+1]//2)
        else:
            chromosome[i] = '-' + str(cycle[2*i]//2)

    return chromosome


def black_edges(genome):
    edges = []

    check = any(isinstance(x, list) for x in genome)
    if check == False:
        l = len(genome)
        nodes = [0]*2*l
        for i in range(l):
            try:
                block = genome[i]
                num = int(block[1:])
                if block[0] == '+':
                    nodes[2*i] = num*2 - 1
                    nodes[2*i+1] = num*2
                elif block[0] == '-':
                    nodes[2*i] = 2*num
                    nodes[2*i+1] = 2*num - 1
            except ValueError:
                break
        for i in range(0, len(nodes), 2):
            edges.append([nodes[i], nodes[i + 1]])

    elif check == True:
        for chromosome in genome:
            l = len(chromosome)
            nodes = [0] * 2 * l
            for i in range(l):
                try:
                    block = chromosome[i]
                    num = int(block[1:])
                    if block[0] == '+':
                        nodes[2 * i] = num * 2 - 1
                        nodes[2 * i + 1] = num * 2
                    elif block[0] == '-':
                        nodes[2 * i] = 2 * num
                        nodes[2 * i + 1] = 2 * num - 1
                except ValueError:
                    break

            for i in range(0, len(nodes), 2):
                edges.append([nodes[i], nodes[i + 1]])

    return edges



def colored_edges(genome): # graphical edges connected synteny blocks
    edges = []
    check = any(isinstance(x, list) for x in genome)

    if check == True:
        for chrom in genome:
            nodes = chrom_to_cycle(chrom)
            l = len(nodes)
            for i in range(1, l-1, 2):
                edges.append([nodes[i], nodes[i+1]])
            edges.append([nodes[-1], nodes[0]])

    elif check == False:
        nodes = chrom_to_cycle(genome)
        l = len(nodes)
        for i in range(1, l - 1, 2):
            edges.append([nodes[i], nodes[i + 1]])
        edges.append([nodes[-1], nodes[0]])

    return edges



from copy import deepcopy

def two_break_graph_to_genome(graph):
    genome = []
    avail_edges = deepcopy(graph)
    cycles = []
    cycle = []

    while len(avail_edges) > 0:
        next = avail_edges[0]
        avail_edges.remove(next)
        cycle.append(next)
        start_node = cycle[0][0]
        end_node = cycle[-1][1]
        while start_node != end_node:
            for edge in avail_edges:
                if end_node in edge:
                    cycle.append(edge)
                    avail_edges.remove(edge)
                    if edge[0] == end_node:
                        end_node = edge[1]
                    elif edge[1] == end_node:
                        end_node = edge[0]
        cycles.append(cycle)
        cycle = []
    for cycle in cycles:
        cycle = cycle[::-1]
        nodes = []
        for i in range(len(cycle)):
            if abs(cycle[i][0] - cycle[i][1]) == 1:
                if cycle[i][0] not in nodes and cycle[i][1] not in nodes:
                    nodes.append(cycle[i][1])
                    nodes.append(cycle[i][0])

        chromosome = cycle_to_chrom(nodes)
        genome.append(chromosome)

    return genome


def cycle_form(graph):
    avail_edges = deepcopy(graph)
    cycles = []
    cycle = []

    while len(avail_edges) > 0:
        next = list(avail_edges[0])
        avail_edges.remove(next)
        cycle.append(list(next))
        start_node = cycle[0][0]
        end_node = cycle[-1][1]
        while start_node != end_node:
            for edge in avail_edges:
                if end_node in edge:
                    cycle.append(edge)
                    avail_edges.remove(edge)
                    if edge[0] == end_node:
                        end_node = edge[1]
                    elif edge[1] == end_node:
                        end_node = edge[0]
        cycles.append(cycle)
        cycle = []

    return cycles



def break_sort(q, p):
    q_form = '(' + ' '.join(b for b in q) + ')'
    print('Original Permutation: ' + q_form)
    print('Genomic Rearrangements:')
    blue_edge = colored_edges(p)
    red_edge = colored_edges(q)
    black_edge = black_edges(p)
    graph = red_edge + blue_edge

    check = any(isinstance(x, list) for x in p)
    if check == True:
        blocks = sum(len(c) for c in p)
    elif check == False:
        blocks = len(p)
    cycles = cycle_form(graph)

    while len(cycles) != blocks:
        nodes = []
        for cycle in cycles:
            if len(cycle) > 2:
                for i in range(len(cycle)-1):
                    if cycle[i] in blue_edge:
                        if cycle[i-1] not in blue_edge and cycle[i+1] not in blue_edge:
                            edge1 = cycle[i-1]
                            edge2 = cycle[i+1]
                            nodes = edge1 + edge2
                            for node in cycle[i]:
                                nodes.remove(node)
                            red_edge.remove(edge1)
                            red_edge.remove(edge2)
                            red_edge.append(cycle[i])
                            red_edge.append(nodes)
                            break

            if len(nodes) > 0:
                break
            else:
                continue

        q = two_break_graph_to_genome(red_edge + black_edge)
        graph = blue_edge + red_edge
        cycles = cycle_form(graph)
        q_form = ''
        for c in q:
            q_form += '(' + ' '.join(b for b in c) + ')'
        print(q_form)

'''
Example
'''



break_sort(['+1', '-2', '-3', '+4'], ['+1', '+2', '-4', '-3'])


# Output: Original Permutation: (+1 -2 -3 +4)
# Genomic Rearrangements:
# (-2 +4 +3)(-1)
# (+3 -1 -2 +4)
# (-2 -1 +3 +4)