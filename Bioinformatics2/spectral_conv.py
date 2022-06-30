'''
Spectral Convolution Solution
'''

'''
A spectral convolution represents all positive mass differences of a given mass spectrum.
The formation of a spectral convolution can be used to account for missing masses in a spectrum
as well as masses that do not conform to a standard set of options. As is the case with
cyclic non-ribosomal peptides. This function returns the m most most commonly occurring masses
in the spectral convolution by sorting a dictionary by values that represent mass multiplicity
and creating separate lists of the masses and their multiplicity in the same index position.
All multiplicities following the m-th mass are compared to ensure ties are kept and the most common
masses are returned as a list.
'''


def spectral_conv(spectrum, m):
    l = len(spectrum)
    conv = {}
    for i in range(l):
        for mass in spectrum:
            diff = int(spectrum[i]) - int(mass)
            if diff > 56:
                if diff < 201:
                    if diff in conv:
                        conv[diff] += 1
                    else:
                        conv.setdefault(diff, 0)
                        conv[diff] += 1
    sorted_conv = dict(sorted(conv.items(), key=lambda x: x[1], reverse=True))
    mass_list = list(sorted_conv.keys())
    multiplicity = list(sorted_conv.values())
    k = len(mass_list)
    if k < int(m):
        return mass_list
    cutoff = multiplicity[m-1]
    for i in range(m, k):
        if multiplicity[i] < cutoff:
            return 'The top ' + str(m) + ' most common masses with ties are: ' + str(mass_list[:i])
    return 'The top ' + str(m) + ' most common masses with ties are: ' + str(mass_list)




'''
Example
'''


print(spectral_conv([0, 57, 118, 179, 236, 240, 301], 3))


# Output: 'The top 3 most common masses with ties are: [61, 122, 57, 118, 179, 183]'
