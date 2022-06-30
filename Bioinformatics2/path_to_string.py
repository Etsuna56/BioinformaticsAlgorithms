'''
String Reconstruction from Path Solution
'''

'''
This function assembled a DNA strand from its complete list of kmers. For this specific function,
the kmers must be given in genomic order. To assemble to genome, the last nucleotide of every
kmer is added to an increasing string that begins with the first kmer in the sequence.
'''


def path_to_str(pattern):
    string = ''
    string += str(pattern[0])
    l = len(pattern)
    for i in range(1, l):
        string += str(pattern[i][-1])
    return string



'''
Example
'''

print(path_to_str(['ACGA', 'ATAC', 'ATCG', 'CATC', 'CGAG', 'CGAT', 'GAGA', 'GATA', 'TACG', 'TCGA']))
# Output: 'ACGACGCGTAAGA'

