'''
Theoretical Spectrum Solution
'''

'''
Utilizing the integer masses of the peptides encoded by ribosomes in the central dogma of biology,
these functions both construct the theoretical spectrum of a given peptide sequence. A linear spectrum
is constructed by finding all possible prefix masses and finding the difference with all possible
suffix masses. For a circular spectrum, the function also adds masses of peptides wrapping around each
linear end of the sequence.
'''



genetic_code = {'AAA': ['K'], 'AAC': ['N'], 'AAG': ['K'], 'AAU': ['N'], 'ACA': ['T'], 'ACC': ['T'], 'ACG': ['T'], 'ACU': ['T'], 'AGA': ['R'], 'AGC': ['S'], 'AGG': ['R'], 'AGU': ['S'], 'AUA': ['I'], 'AUC': ['I'], 'AUG': ['M'], 'AUU': ['I'], 'CAA': ['Q'], 'CAC': ['H'], 'CAG': ['Q'], 'CAU': ['H'], 'CCA': ['P'], 'CCC': ['P'], 'CCG': ['P'], 'CCU': ['P'], 'CGA': ['R'], 'CGC': ['R'], 'CGG': ['R'], 'CGU': ['R'], 'CUA': ['L'], 'CUC': ['L'], 'CUG': ['L'], 'CUU': ['L'], 'GAA': ['E'], 'GAC': ['D'], 'GAG': ['E'], 'GAU': ['D'], 'GCA': ['A'], 'GCC': ['A'], 'GCG': ['A'], 'GCU': ['A'], 'GGA': ['G'], 'GGC': ['G'], 'GGG': ['G'], 'GGU': ['G'], 'GUA': ['V'], 'GUC': ['V'], 'GUG': ['V'], 'GUU': ['V'], 'UAA': ['0'], 'UAC': ['Y'], 'UAG': ['0'], 'UAU': ['Y'], 'UCA': ['S'], 'UCC': ['S'], 'UCG': ['S'], 'UCU': ['S'], 'UGA': ['0'], 'UGC': ['C'], 'UGG': ['W'], 'UGU': ['C'], 'UUA': ['L'], 'UUC': ['F'], 'UUG': ['L'], 'UUU': ['F']}

mass_code = {'G': 57, 'A': 71, 'S': 87, 'P': 97, 'V': 99, 'T': 101, 'C': 103, 'I': 113, 'L': 113, 'N': 114, 'D': 115, 'K': 128, 'Q': 128, 'E': 129, 'M': 131, 'H': 137, 'F': 147, 'R': 156, 'Y': 163, 'W': 186}



def linear_spec(peptide):
    mass_dict = mass_code
    l = len(peptide)
    pref_mass = [0]
    for i in range(l):
        for pep in mass_dict:
            if pep == peptide[i]:
                pref_mass.append(pref_mass[i] + mass_dict[pep])
    lin_spec = [0]
    for i in range(l):
        for j in range(i+1, l+1):
            lin_spec.append(int(pref_mass[j] - pref_mass[i]))
    return 'Theoretical linear mass spectrum: ' + str(sorted(lin_spec))



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
    return 'Theoretical cyclic mass spectrum: ' + str(sorted(circ_spec))



'''
Example
'''


print(linear_spec('NQEL'))

# Output: 'Theoretical linear mass spectrum: [0, 113, 114, 128, 129, 242, 242, 257, 370, 371, 484]'


print(cyclic_spec('NQEL'))

# Output: 'Theoretical cyclic mass spectrum: [0, 113, 114, 128, 129, 227, 242, 242, 257, 355, 356, 370, 371, 484]'