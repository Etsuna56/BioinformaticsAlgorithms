'''
Greedy Genome Sorting Solution
'''

'''
The following algorithm outputs the sequence of sorting operations required to transform a given permutation
into an identity permutation with an equal number of synteny blocks. This is done by iteratively checking whether
every block is correctly positioned and performing a reversal in the case that a block is out of position. This
is continued until every element in the permutation is sorted, at which point the total number of operations required
is returned.
'''




def greed_sorting(perm):
    p = len(perm)
    identity_perm = [('+' + str(i)) for i in range(1, p + 1)]
    count = 0
    print('Sorting Operations:')

    while perm != identity_perm:
        for sb in range(1, p + 1):
            if perm[sb - 1] != identity_perm[sb - 1]:
                gen_rev = []
                k = perm.index('+' + str(sb)) if ('+' + str(sb)) in perm else perm.index('-' + str(sb))
                gen_rev = list(reversed(perm[sb - 1:k + 1]))
                for i in range(len(gen_rev)):
                    block = gen_rev[i]
                    if block[0] == '-':
                        block = block.replace('-', '+')
                        gen_rev[i] = block
                    else:
                        block = block.replace('+', '-')
                        gen_rev[i] = block
                perm[sb - 1:k + 1] = gen_rev
                print(*perm, sep=' ')
                count += 1
                if perm[sb - 1] != identity_perm[sb - 1]:
                    perm[sb - 1] = perm[sb - 1].replace('-', '+')
                    print(*perm, sep=' ')
                    count += 1

    print('Number of Sorting Operations: ' + str(count))


'''
Example
'''



greed_sorting(['+6','-4', '+5', '-3','+8', '+2', '+7','-1'])

# Output: Sorting Operations:
# +1 -7 -2 -8 +3 -5 +4 -6
# +1 +2 +7 -8 +3 -5 +4 -6
# +1 +2 -3 +8 -7 -5 +4 -6
# +1 +2 +3 +8 -7 -5 +4 -6
# +1 +2 +3 -4 +5 +7 -8 -6
# +1 +2 +3 +4 +5 +7 -8 -6
# +1 +2 +3 +4 +5 +6 +8 -7
# +1 +2 +3 +4 +5 +6 +7 -8
# +1 +2 +3 +4 +5 +6 +7 +8
# Number of Sorting Operations: 9