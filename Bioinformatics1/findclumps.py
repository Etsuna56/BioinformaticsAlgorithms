# ''''''
# Clump Finder Solution
# ''''''



def frequencytable(text, k):
    freq_map = {}
    n = len(text)
    for i in range(n - k + 1):
        pattern = text[i:i + k]
        if pattern in freq_map:
            freq_map[pattern] += 1
        if pattern not in freq_map:
            freq_map[pattern] = 1
    return freq_map




def findclumps(text, k, l, t): # Returns kmers with at least t repetitions within substrings of length l
    # text = entire sequence
    # k = length of kmer
    # l = length of window within text
    # t = minimum repetitions of kmer within l
    patterns = []
    count = 0
    n = len(text)
    for i in range(n - l + 1):
        window = text[i:i + l]
        freq_map = frequencytable(window, k)
        for pattern in freq_map:
            if freq_map[pattern] >= t:
                if pattern not in patterns:
                    patterns.append(pattern)
                    count += 1
    return patterns, count

# ''''''
# Example
# ''''''

print(findclumps('ACGCGCGATCGCGATCTACGATCGCGCG', 2, 8, 3))
# output: (['CG'], 1)
