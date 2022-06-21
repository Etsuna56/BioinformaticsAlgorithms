# ''''''
# Hamming Distance Solution
# ''''''



def hammingdistance(pattern, string): # Returns a count for nucleotide differences between strings
    dist = 0
    for i in range(len(pattern)):
        if pattern[i] != string[i]:
            dist += 1
    return dist


# ''''''
# Example
# ''''''

print(hammingdistance('ACAACGCCGCGGGCGTTG', 'ACGCATCGACTAGTTCGAT'))
# output: 13