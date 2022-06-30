'''
Eulerian Path Solution
'''

'''
A Eulerian path in a deBruijn graph represents a path through the graph that crosses every edge
exactly once. In order for a graph to have a Eulerian path, at most two nodes can be unbalanced.
Of these two nodes, the node with more outward edges will be used as the starting position, and
the node with more inward edges will be the final node. This algorithm searches for a Eulerian
path in a deBruijn graph by first determining which nodes are unbalanced by counting their outward
edges and comparing this to all other nodes. Once a start position has been determined, a depth
first search is performed until a path is found that exhausts all edges, and the final path through
the graph is returned with directed edges connecting nodes in sequence.
'''



from random import randint
from copy import deepcopy


def create_adj_list(graph):
    adj_list = {} # key: node, value: nodes connected by outward edges
    circuit_max = 0
    for line in graph:
        node = line.strip('\n') # form: 1: 2 3
        node = node.replace(':', '') # form: 1 2 3
        node = node.replace(' ', ',') # form: 1,2,3
        node = node.replace(',', ' ', 1) # form: 1 2,3
        node = node.split(' ') # form: 1 \n 2,3
        adj_list.setdefault(node[0], [])
        for number in node[1].split(','):
            adj_list[node[0]].append(number) # adj_list['1] ->  '1': ['2', '3']
            circuit_max += 1
    return adj_list, circuit_max



def find_start(red_adj_list):
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
    return red_adj_list, start_node


def eulerian_path(graph):
    # If inporting from file:
    # adj_list, circuit_max = create_adj_list(graph)
    # If using dictionary form:
    adj_list = graph
    circuit_max = sum(len(i) for i in list(adj_list.values()))
    red_adj_list = {}
    red_adj_list = deepcopy(adj_list)
    red_adj_list, start_node = find_start(red_adj_list)
    start = start_node
    curr_vert = start_node
    stack = []
    circuit = []  # eulerian path build as nodes run out of edges
    while len(circuit) != circuit_max:  # continues increasing length of path until every edge has been crossed
        if len(red_adj_list[curr_vert]) > 0:  # checks if outward edges are available
            stack.append(curr_vert)
            pick = randint(0, len(red_adj_list[curr_vert]) - 1)  # picks random remaining outward edge
            temp = deepcopy(curr_vert)  # temporary copy of current vertex for indexing
            curr_vert = red_adj_list[temp][pick]  # random node to continue current path
            red_adj_list[temp].remove(curr_vert)  # remove edge connecting current and previous nodes
        else:
            circuit.append(curr_vert)  # adds node to final path once it has no outward edges left
            curr_vert = stack[len(stack) - 1]
            stack.pop() # removes from list of nodes with available edges
    path = ''
    path += start
    for vert in reversed(circuit):
        path += ('->' + vert)
    return path




'''
Example
'''


print(eulerian_path({'TA': ['AA', 'AG'], 'AA': ['AT'], 'AT': ['TG'], 'TG': ['GT'], 'GT': ['TT'], 'TT': ['TA']}))


# Output: TA->AA->AT->TG->GT->TT->TA->AG
