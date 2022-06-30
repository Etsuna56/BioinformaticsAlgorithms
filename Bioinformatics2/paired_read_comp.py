'''
Read Pair Composition Solution
'''

'''
During genome sequencing, the length of reads if an important metric to optimize the assembly
process. Read pairs are one method that biologists employ to obtain more information from a
collection of DNA reads. This function forms all read pairs of a given DNA string by concatenating
two kmers with a given distance between them. This is repeated for all possible kmers with d distance
between the two kmers. The final set of paired reads is returned in lexicographical order because
this is more reasonable for genome sequencing since the order of the original genome is unknown.
'''



def pairedcomp(text, k, d):
    paired_reads = []
    for i in range(len(text) - 2*k - d + 1):
        read = ''
        if len(text[i+k+d:i+d+(k*2)]) == k:
            read += text[i:i+k] + '|' + text[i+k+d:i+d+(k*2)]
        else:
            break
        paired_reads.append(read)
    return sorted(paired_reads)


'''
Example
'''


print(pairedcomp('AGCCTAACGCATC', 3, 2))

# Output: ['AAC|ATC', 'AGC|AAC', 'CCT|CGC', 'CTA|GCA', 'GCC|ACG', 'TAA|CAT']