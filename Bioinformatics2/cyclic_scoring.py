'''
Cyclic Peptide Scoring Solution
'''

'''
These functions compute a score for a given peptide sequence and an experimental spectrum
by comparing each element of the experimental spectrum to each element in a theoretical
spectrum for the peptide sequence and adding a point to the score for every element they share.
'''


mass_code = {'G': 57, 'A': 71, 'S': 87, 'P': 97, 'V': 99, 'T': 101, 'C': 103, 'I': 113, 'L': 113, 'N': 114, 'D': 115, 'K': 128, 'Q': 128, 'E': 129, 'M': 131, 'H': 137, 'F': 147, 'R': 156, 'Y': 163, 'W': 186}

spec_masses = {57: 'G', 71: 'A', 87: 'S', 97: 'P', 99: 'V', 101: 'T', 103: 'C', 113: 'L', 114: 'N', 115: 'D', 128: 'Q', 129: 'E', 131: 'M', 137: 'H', 147: 'F', 156: 'R', 163: 'Y', 186: 'W'}



def cyclic_spec(peptide):
    mass_dict = mass_code
    l = len(peptide)
    pref_mass = [0]
    for i in range(l):
        for pep in mass_dict:
            if pep == peptide[i]:
                pref_mass.append(pref_mass[i] + mass_dict[pep])
    circ_spec = [0]
    peptide_mass = pref_mass[l]
    for i in range(l):
        for j in range(i+1, l+1):
            circ_spec.append(pref_mass[j] - pref_mass[i])
            if i > 0 and j < l:
                circ_spec.append(peptide_mass - (pref_mass[j] - pref_mass[i]))
    return sorted(circ_spec)


def cyclopep_score(peptide, spectrum):
    theory_spec = cyclic_spec(peptide)
    score = 0
    x = 0
    for i in range(len(theory_spec)):
        for j in range(x, len(spectrum)):
            if spectrum[j] == str(theory_spec[i]):
                score += 1
                x = j + 1
                break
    return 'The spectra have ' + str(score) + ' masses in common.'


'''
Example
'''


print(cyclopep_score('NQEL', ['0', '57', '57', '113', '114', '128', '129', '242', '242', '257', '370', '371', '484']))


# Output: 'The spectra have 11 masses in common'
