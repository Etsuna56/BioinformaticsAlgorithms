# ''''''
# Skew Solution
# ''''''



def skew(text): # Returns a positional skew which represents the variation in occurrence of C relative to G
    skew_count = 0
    skew_list = [0]
    for i in range(len(text)):
        if text[i] == 'G':
            skew_count += 1
        elif text[i] == 'C':
            skew_count -= 1
        else:
            skew_count = skew_count
        skew_list.append(skew_count)
    return skew_list


# ''''''
# Example
# ''''''

print(skew('ACAACGCCGCGGGCGTTG'))
# output: [0, 0, -1, -1, -1, -2, -1, -2, -3, -2, -3, -2, -1, 0, -1, 0, 0, 0, 1]
