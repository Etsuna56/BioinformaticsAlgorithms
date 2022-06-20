# ''''''
# Frequency Table Solution
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



# ''''''
# Example
# ''''''

print(frequencytable('AGCTATCGCGCAATC', 3))
# output: {'AGC': 1, 'GCT': 1, 'CTA': 1, 'TAT': 1, 'ATC': 2, 'TCG': 1, 'CGC': 2, 'GCG': 1, 'GCA': 1, 'CAA': 1, 'AAT': 1}
