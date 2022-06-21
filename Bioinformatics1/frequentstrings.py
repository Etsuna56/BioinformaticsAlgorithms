# ''''''
# Pattern Frequency Solution
# ''''''


def frequent_strings(text, k):
    frequent_patterns = {}
    for i in range(len(text) - k + 1):
        pattern = text[i:i + k]
        if pattern in frequent_patterns:
            frequent_patterns[pattern] += 1
        if pattern not in frequent_patterns:
            frequent_patterns[pattern] = 1
    return frequent_patterns


# ''''''
# Example
# ''''''

print(frequent_strings('CGGAGGACTTCA', 3))
# output: {'CGG': 1, 'GGA': 2, 'GAG': 1, 'AGG': 1, 'GAC': 1, 'ACT': 1, 'CTT': 1, 'TTC': 1, 'TCA': 1}
