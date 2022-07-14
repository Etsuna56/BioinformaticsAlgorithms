'''
Longest Path in DAG Solution
'''

'''
The following algorithm solves the longest path in a directed acyclic graph, G, by instead finding the shortest
path in graph G' where all edge weights are converted to their negative value. The scores of all nodes except
the source are first initialized as infinity. A topological ordering is determined by ensuring all nodes with
paths connecting to subsequent nodes are visited before moving on. This order ensures that every incoming edge
has been considered when determining the shortest path to a node. Once all shortest paths have been determined,
the longest path to any node in the original graph is simply the negative of the shortest path. The nodes visited
in order to form this path are also recorded by determining which edges followed in the graph yielded the resulting
score.
'''


import sys


def add_edge(u, v, w):
    global adj

    adj[u].append([v, w])

def form_path(v):
    global path, prev_node

    path.append(v)
    if len(prev_node[v]) > 0:
        form_path(prev_node[v][0])

def longest_path(v):
    global visited, adj, top_order, prev_node, path, nodes

    for node in adj[v]:
        try:
            if visited[node[0]] == 0:
                longest_path(node[0])
        except IndexError:
            print(node[0], 'caused index error')
            break

    if v not in top_order:
        top_order.append(v)

def final_path(s):
    global visited, top_order, adj, n, path, nodes, prev_node, e

    dist = [sys.maxsize]*n
    dist[s] = 0

    for i in range(n):
        if visited[i] == 0: # forms topological order
            longest_path(i)

    while len(top_order) > 0:
        u = top_order[-1] # follows topological order and finds longest path to every node in order
        del top_order[-1]
        if dist[u] != sys.maxsize:
            for v in adj[u]:
                if dist[v[0]] > dist[u] + v[1] * -1:
                    dist[v[0]] = dist[u] + v[1] * -1
                    if u not in prev_node[v[0]]:
                        prev_node[v[0]] = [u] # previous node in longest path to node v

    for i in range(n):
        if dist[i] == min(dist):
            if i == e:
                form_path(i) # forms longest to end node

    print('Longest path score: ' + str(min(dist) * -1))
    ans = ''
    for node in path[::-1]:
        if node != path[0]:
            ans += str(node) + '->'
        else:
            ans += str(node)

    print('Longest Path: ' + ans)


'''
Example
'''


if __name__ == '__main__':
    x = [0,0,0,1,1,2,2,2,3,3,4,5]
    y = [1,2,3,2,5,4,5,6,4,5,6,6]
    z = [5,6,5,2,4,4,3,5,6,8,2,1]

    nodes = set()
    for i in x:
        nodes.add(i)
    for i in y:
        nodes.add(i)

    n = len(nodes) + 1

    visited = [0 for i in range(n)]
    top_order = []
    adj = [[] for i in range(n)]
    prev_node = {}
    path = []

    for i in range(n):
        prev_node.setdefault(i, [])

    for u,v,w in zip(x,y,z):
        add_edge(u, v, w)

    s = 0
    e = max(nodes)
    final_path(s)


# Output: Longest path score: 14
# Longest Path: 0->3->5->6