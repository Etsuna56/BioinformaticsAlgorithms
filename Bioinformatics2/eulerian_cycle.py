'''
Eulerian Cycle Solution
'''

'''
A Eulerian cycle in a deBruijn graph represents a path through the graph that path through the graph
that crosses every edge while ending at the same start point. In order for a graph to have a
Eulerian circuit, all nodes must be balances and the graph must be well connected, meaning it must
be possible to reach any node from any other node. Given a adjacency list in the form of a dictionary,
this algorithm searches for a Eulerian cycle through the graph by performing a depth first search until
every edge has been traversed, which occurs once the final circuit matches the maximum circuit length.
This final circuit is then reversed to correct the order and the path through the graph is returned.
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



def eulerian_cycle(graph):
    # If inporting from file:
    # adj_list, circuit_max = create_adj_list(graph)
    # If using dictionary form:
    adj_list = graph
    circuit_max = sum(len(i) for i in list(adj_list.values()))
    red_adj_list = {}
    red_adj_list = deepcopy(adj_list)
    start = list(red_adj_list.keys())[0] # arbitrary start position
    curr_vert = list(red_adj_list.keys())[0]
    stack = []
    circuit = []
    while len(circuit) != circuit_max: # continues increasing length of path until every edge has been crossed
        if len(red_adj_list[curr_vert]) > 0: # checks if outward edges are available
            stack.append(curr_vert)
            pick = randint(0, len(red_adj_list[curr_vert]) - 1) # picks random outward edge
            temp = deepcopy(curr_vert) # temporary copy of current vertex for indexing
            curr_vert = red_adj_list[temp][pick] # random node to continue current path
            red_adj_list[temp].remove(curr_vert) # remove edge connecting current and previous nodes
        else:
            circuit.append(curr_vert)
            curr_vert = stack[len(stack) - 1]
            stack.pop() # removes node from current path once it has no outward edges left
    path = ''
    path += start
    for vert in reversed(circuit):
        path += ('->' + vert)
    return path


'''
Example
'''


print(eulerian_cycle({'TA': ['AA'], 'AA': ['AT'], 'AT': ['TG'], 'TG': ['GT'], 'GT': ['TT'], 'TT': ['TA']}))


# Output: TA->AA->AT->TG->GT->TT->TA
