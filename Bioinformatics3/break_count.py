'''
Breakpoint Count Solution
'''

'''
Given a permutation of a genome, a breakpoint represents a point in the permutation where the order of elements
breaks from the numerical ordering expected in an identity permutation. The following algorithm computed the number
of breakpoints in a permutation by comparing each i-th element to the i+1 element and determining the difference
between the elements to identify them as a breakpoint or an adjacency.
'''




def break_count(perm):
    perm.insert(0, '0')
    p = len(perm)
    perm.insert(p, '+' + str(p))
    count = 0

    for i in range(p):
        if perm[i][0] == '-':
            if int(perm[i+1][1:]) - int(perm[i][1:]) != -1:
                count += 1
        elif perm[i][0] == '+':
            if int(perm[i+1][1:]) - int(perm[i][1:]) != 1:
                count += 1
        elif perm[i] == '0':
            if perm[i+1] != '+1':
                count += 1

    print('Breakpoint Count: ' + str(count))


'''
Example
'''



break_count(['+6','-4', '+5', '-3','+8', '+2', '+7','-1'])

# Output: Breakpoint Count: 9