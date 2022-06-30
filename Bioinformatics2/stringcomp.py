'''
String Composition Solution
'''

'''
Given a DNA strand as a string, this function returns a list of all possible kmers formed
from the DNA strand. This is performed by sliding a window through the DNA strand and
recording each k length window as a separate kmer.
'''


def str_comp(text, k):
    kmers = []
    l = len(text)
    for i in range(l - k + 1):
        kmer = text[i:i + k]
        kmers.append(kmer)
    kmers = sorted(kmers)
    return kmers



'''
Example
'''

print(str_comp('CATCGATACGAGA', 4))
# Output: ['ACGA', 'ATAC', 'ATCG', 'CATC', 'CGAG', 'CGAT', 'GAGA', 'GATA', 'TACG', 'TCGA']

