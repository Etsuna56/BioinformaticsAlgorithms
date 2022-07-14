'''
2-Break Distance Solution
'''

'''
The 2-break distance between two genomes is the number of genomic rearrangements required to transform one
genome into another. This value is equivalent to the number of synteny blocks in a breakpoint graph between the two
genomes minus the number of alternating cycles within this breakpoint graph before any rearrangements are performed.
These alternating cycles are alternating edges that correspond to the formation of the each genome within the graph.
These cycles are found by forming continuous sequences of edges until the final node of the sequence is equivalent to
the starting sequence, at which point the cycle would repeat. The total number of cycles in the graph will be found
once all available edges in the breakpoint graph have been assigned to cycles.
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


def graph_to_genome(graph): # graphical colored edges split into cycles and converted to chromosomal permutations
    genome = []
    cycles = []
    cycle = []
    l = len(graph)

    for i in range(l-1):
        if abs(graph[i][1] - graph[i+1][0]) == 1:
            if graph[i] not in cycle:
                cycle.append(graph[i])
                cycle.append(graph[i+1])
                if graph[i+1] == graph[-1]:
                    cycles.append(cycle)
            else:
                cycle.append(graph[i+1])
                if graph[i+1] == graph[-1]:
                    cycles.append(cycle)
        else:
            if len(cycle) > 0:
                cycles.append(cycle)
                cycle = [graph[i+1]]

    for cycle in cycles:
        nodes = []
        for i in range(len(cycle)):
            nodes.append(cycle[i-1][1])
            nodes.append(cycle[i][0])

        chromosome = cycle_to_chrom(nodes)
        genome.append(chromosome)

    return genome


from copy import deepcopy

def break_dist(p, q):
    genome = []
    blocks = 0

    check = any(isinstance(x, list) for x in q)
    if check == True:
        for i in range(len(q)):
            genome.append(q[i])
    elif check == False:
        genome.append(q)

    check = any(isinstance(x, list) for x in p)
    if check == True:
        # print('p is multichromal')
        for i in range(len(p)):
            blocks += len(p[i])
            genome.append(p[i])

    elif check == False:
        # print('p is unichromal')
        blocks = len(p)
        genome.append(p)

    graph = colored_edges(genome)
    avail_edges = deepcopy(graph)
    l = len(graph)
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

    print('Breakpoint Distance: ' + str((blocks - len(cycles))))

    return (blocks - len(cycles))


'''
Example
'''



break_dist(['+1', '+2', '+3', '+4', '+5', '+6'], [['+1', '-3', '-6', '-5'], ['+2', '-4']])


# Output: Breakpoint Distance: 3