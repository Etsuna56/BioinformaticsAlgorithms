# ''''''
# Reverse Complement String Solution
# ''''''




def reversecomplement(strand): # Takes a DNA strand and returns its complementary strand
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



# ''''''
# Example
# ''''''

print(reversecomplement('AGCTATCGCGCAATC'))
# output: 'GATTGCGCGATAGCT'
