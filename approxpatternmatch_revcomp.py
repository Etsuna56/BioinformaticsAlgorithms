# ''''''
# Pattern Match with Hamming Distance and Reverse Complement Solution
# ''''''


# This series of functions identifies matches of a pattern within larger strings with at most d mismatches
# while also accounting for the double-stranded nature of DNA by checking for reverse complement patterns within
# the same string.




def reversecomplement(strand):
    complement = ''
    for i in strand:
        if i == 'A':
            complement += 'T'
        elif i == 'T':
            complement += 'A'
        elif i == 'C':
            complement += 'G'
        elif i == 'G':
            complement += 'C'
    return complement[::-1]




def hammingdistance(pattern, string): # Returns a count for nucleotide differences between strings
    dist = 0
    for i in range(len(pattern)):
        if pattern[i] != string[i]:
            dist += 1
    return dist



def freqwords_mismatch_revcomp(text, string, d):
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
    for i in range(len(text) - k + 1):
        pattern = text[i:i + k]
        rev_str = reversecomplement(string)
        if hammingdistance(pattern, rev_str) <= d:
            if pattern in matches:
                matches[pattern] += 1
                count += 1
                match_pos.append(i)
            if pattern not in matches:
                matches[pattern] = 1
                count += 1
                match_pos.append(i)
    matches = sorted(matches.items(), key=lambda x: x[1], reverse=True)
    print('matches: ' + str(matches) + '\n' 'match locations: ' + str(match_pos) + '\n' 'number of matches: ' + str(count))


# ''''''
# Example
# ''''''



freqwords_mismatch_revcomp('GCCACGACGATCTAGATCGACGATCTACGGGCT', 'AGCT', 1)
# matches: [('ATCT', 4), ('AGAT', 2), ('GGCT', 2)]
# match locations: [9, 13, 22, 29, 9, 13, 22, 29]
# number of matches: 8