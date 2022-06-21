# ''''''
# String Pattern Match Solution
# ''''''


import re



def patternmatch(text, pattern): # Returns starting positions for a substring within larger strings as a list
    matches = []
    for match in re.finditer('(?=' + pattern + ')', text):
        matches.append(match.start())
    return matches


# ''''''
# Example
# ''''''

print(patternmatch('AGCTATCGCGCAATC', 'AT'))
# output: [4, 12]
