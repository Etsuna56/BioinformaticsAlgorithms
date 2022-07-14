'''
Manhattan Tourist with Diagonal Edges Solution
'''

'''
The Manhattan Tourist problem with Diagonal Edges represents the problem of finding the longest path 
through a two dimensional weighted grid. In this grid, one must move south, east, or southeast at every node.
To solve this problem, the following algorithm builds a scoring matrix that records the longest
path to every node in the grid. This is done by starting at the source node and understanding 
that every subsequent node's longest path is the maximum of the possible paths leading to this node. 
In the case of the two dimensional Manhattan grid, each node has a path from the north, west, and northwest.
So each nodes longest path of the maximum of the weighted edges leading to the node plus the score of the node
from which the edge left.
'''


import numpy as np


def manhat_tourist_diag(n, m, down, right, diag):

    s = np.matrix(np.arange((n+1)*(m+1)).reshape((n+1, m+1)))
    s[0, 0] = 0

    for i in range(1, n + 1):
        s[i, 0] = s[i-1, 0] + down[i-1, 0]
    for j in range(1, m + 1):
        s[0, j] = s[0, j-1] + right[0, j-1]
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            s[i, j] = max([s[i-1, j] + down[i-1, j], s[i, j-1] + right[i, j-1], s[i-1, j-1] + diag[i-1, j-1]])

    print('Score Matrix: \n' + str(s))
    return 'Longest Path to Sink: ' + str(s[n, m])


'''
Example
'''

right = np.matrix([[3,2,4,0], [3,2,4,2], [0,7,3,4], [3,3,0,2], [1,3,2,2]])

down = np.matrix([[1,0,2,4,3], [4,6,5,2,1], [4,4,5,2,1], [5,6,8,5,3]])

diag = np.matrix([[5,0,2,1], [8,4,3,0], [10,8,9,5], [5,6,4,7]])

print(manhat_tourist_diag(4, 4, down, right, diag))


# Output: Score Matrix:
# [[ 0  3  5  9  9]
#  [ 1  5  7 13 15]
#  [ 5 11 18 21 25]
#  [ 9 15 23 27 29]
#  [14 21 31 33 35]]
# Longest Path to Sink: 35
