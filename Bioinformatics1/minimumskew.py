# ''''''
# Minimum Skew Solution
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


def min_skew(text): # Returns minimum skew positions, which represent origins of replication in bacteria
    skew_list = skew(text)
    minimum = min(skew_list)
    min_skew_list = []
    for i in range(len(text)):
        if skew_list[i] == minimum:
            min_skew_list.append(i)
    min_skew_list = str(min_skew_list)
    min_skew_list = min_skew_list.replace(',', '')
    min_skew_list = min_skew_list.replace('[', '')
    min_skew_list = min_skew_list.replace(']', '')
    return min_skew_list


# ''''''
# Example
# ''''''

print(min_skew('ACAACGCCGCGGGCGTTG'))
# output: '8 10'
