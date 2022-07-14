'''
Manhattan Tourist Solution
'''

'''
The Manhattan Tourist problem represents the problem of finding the longest path through a two
dimensional weighted grid. In this grid, one must move either south or east at every node.
To solve this problem, the following algorithm builds a scoring matrix that records the longest
path to every node in the grid. This is done by starting at the source node and understanding 
that every subsequent node's longest path is the maximum of the possible paths leading to this node. 
In the case of the two dimensional Manhattan grid, each node has a path from the north and from the west.
'''


import numpy as np


def manhat_tourist(n, m, down, right):

    s = np.matrix(np.arange((n+1)*(m+1)).reshape((n+1, m+1)))
    s[0, 0] = 0

    for i in range(1, n + 1):
        s[i, 0] = s[i-1, 0] + down[i-1, 0]
    for j in range(1, m + 1):
        s[0, j] = s[0, j-1] + right[0, j-1]
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            s[i, j] = max([s[i-1, j] + down[i-1, j], s[i, j-1] + right[i, j-1]])

    print('Score Matrix: ' + str(s))
    return 'Longest Path to Sink: ' + str(s[n, m])




'''
Example
'''

'''
down = [[0 2 2 1 0 2 3 4 2 1 2 2 0 4 4 2]
 [3 1 0 3 2 2 1 2 1 2 3 3 3 3 4 1]
 [3 0 1 4 1 4 0 1 3 4 0 3 2 1 1 4]
 [4 2 0 3 1 4 1 0 3 3 2 0 0 1 1 3]
 [1 0 1 1 4 1 1 0 4 4 3 2 1 3 0 2]
 [0 3 3 0 3 4 3 1 4 2 1 3 2 2 1 2]
 [3 4 1 4 2 1 0 2 2 3 0 1 0 2 0 2]
 [2 4 3 2 4 2 2 2 4 2 4 4 3 0 4 3]
 [0 1 3 2 1 3 2 4 0 2 2 1 1 1 1 2]
 [1 4 4 4 2 0 0 1 3 3 2 0 1 1 3 4]
 [0 1 1 4 0 0 2 3 1 4 3 0 1 4 0 4]]
'''

'''
right = [[1 1 0 0 4 2 3 2 4 2 4 4 4 4 1]
 [0 3 0 0 0 2 2 3 0 1 1 2 3 3 1]
 [2 4 2 2 4 2 4 3 4 3 1 3 4 3 4]
 [4 0 2 2 4 1 2 0 1 4 3 1 3 4 0]
 [2 2 2 0 0 1 1 2 2 1 1 1 0 0 4]
 [1 1 4 1 3 0 2 4 3 2 1 0 0 3 0]
 [1 3 1 0 0 4 4 1 3 3 0 1 2 4 0]
 [2 0 4 0 1 2 3 2 2 0 2 2 1 3 2]
 [1 2 3 2 4 1 2 4 4 4 4 2 4 4 3]
 [2 1 3 4 4 4 4 2 1 1 3 4 4 2 2]
 [0 2 1 1 2 1 4 0 3 1 3 0 2 4 3]
 [3 1 3 0 4 2 3 2 3 0 2 2 1 0 3]]
'''

'''
n = 11
m = 15
'''


print(manhat_tourist(n, m, down, right))


# Output: Score Matrix: [[ 0  1  2  2  2  6  8 11 13 17 19 23 27 31 35 36]
#  [ 0  3  6  6  6  8 11 15 18 18 21 25 27 35 39 40]
#  [ 3  5  9 11 13 17 19 23 26 30 33 34 37 41 44 48]
#  [ 6 10 10 15 17 21 22 24 29 34 38 41 42 45 49 52]
#  [10 12 14 18 18 25 26 27 32 37 40 41 42 46 50 55]
#  [11 12 15 19 22 26 27 29 36 41 43 44 44 49 52 57]
#  [11 15 18 19 25 30 34 38 40 43 46 47 48 51 55 59]
#  [14 19 19 23 27 31 34 40 42 46 46 48 50 53 56 61]
#  [16 23 25 28 31 35 36 42 46 50 54 58 60 64 68 71]
#  [16 24 28 31 35 39 43 47 49 52 56 59 63 67 69 73]
#  [17 28 32 35 37 39 43 48 52 55 58 61 64 68 72 77]
#  [17 29 33 39 39 43 45 51 53 59 61 63 65 72 72 81]]
# Longest Path to Sink: 81
