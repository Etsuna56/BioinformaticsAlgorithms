'''
K-Universal Binary String Solution
'''

'''
These functions first generate all possible binary strings of k^2 + (k - 1) length without their
identifying prefix. A k length window is then slid through these strings to create all possible
kmers within the string. These kmers are then added to a set, which is compared to the maximum
number of kmers to determine if the string was k-universal.
'''


def gen_strs(n):
    return [bin(i)[2:].zfill(n) for i in range(2**n)]


def universal_str_check(k):
    bin_strs = gen_strs(k**2 + k - 1)
    for string in bin_strs:
        kmers = set()
        for i in range(len(string) - k + 1):
            kmers.add(string[i:i + k])
        if len(kmers) == len(string) - k + 1:
            return string



'''
Example
'''


print(universal_str_check(4))

# Output: '0000100110101111000'

