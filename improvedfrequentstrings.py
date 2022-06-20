# ''''''
# Improved Frequent Strings Solution
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



def freq_str_2(text, k):
    freq_map = frequencytable(text, k)
    freq_map = sorted(freq_map.items(), key=lambda x: x[1], reverse=True)
    return freq_map



# ''''''
# Example
# ''''''

print(freq_str_2('AGCTATCGCGCAATC', 3))
# output: [('ATC', 2), ('CGC', 2), ('AGC', 1), ('GCT', 1), ('CTA', 1), ('TAT', 1), ('TCG', 1), ('GCG', 1), ('GCA', 1), ('CAA', 1), ('AAT', 1)]