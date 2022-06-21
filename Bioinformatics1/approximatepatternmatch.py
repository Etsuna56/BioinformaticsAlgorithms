# ''''''
# Pattern Match with Hamming Distance Solution
# ''''''


# The approximate pattern count function searches through an input string and compares each pattern
# of length k to a given substring and returns patterns with start positions that match the substring
# with at most d nucleotide mismatches.



def hammingdistance(pattern, string): # Returns a count for nucleotide differences between strings
    dist = 0
    for i in range(len(pattern)):
        if pattern[i] != string[i]:
            dist += 1
    return dist



def approxpatterncount(text, string, d):
    k = len(string)
    matches = {}
    match_pos = []
    count = 0
    for i in range(len(text) - k + 1):
        pattern = text[i:i + k]
        if hammingdistance(pattern, string) <= d:
            match_pos.append(i)
            count += 1
            if pattern in matches:
                matches[pattern] += 1
            if pattern not in matches:
                matches[pattern] = 1
    matches = sorted(matches.items(), key=lambda x: x[1], reverse=True)
    print('matches: ' + str(matches) + '\n' 'match locations: ' + str(match_pos) + '\n' 'number of matches: ' + str(count))


# ''''''
# Example
# ''''''

approxpatterncount('ACAACGCCGCGGGCGTTG', 'CGC', 1)
# output: matches: [('CGC', 2), ('CGG', 1), ('GGC', 1), ('CGT', 1)]
# match locations: [4, 7, 9, 11, 13]
# number of matches: 5
