'''
Branch and Bound Cyclopeptide Sequencing Solution
'''

'''
This algorithm is a branch and bound approach to sequencing a peptide from its
mass spectrum. It functions by building all possible candidates from the integer masses
of the 20 ribosomal amino acids one amino acid at a time. After each extension of the
candidate amino acids, they are checked to determine if their theoretical spectrum are
consistent with the experimental spectrum provided. Once a peptide has been built that
matches the largest mass in and is consistent with the original experimental spectrum, it
is added to a final list of peptides, which are returned once all candidate peptides have
been exhausted.
'''



genetic_code = {'AAA': ['K'], 'AAC': ['N'], 'AAG': ['K'], 'AAU': ['N'], 'ACA': ['T'], 'ACC': ['T'], 'ACG': ['T'], 'ACU': ['T'], 'AGA': ['R'], 'AGC': ['S'], 'AGG': ['R'], 'AGU': ['S'], 'AUA': ['I'], 'AUC': ['I'], 'AUG': ['M'], 'AUU': ['I'], 'CAA': ['Q'], 'CAC': ['H'], 'CAG': ['Q'], 'CAU': ['H'], 'CCA': ['P'], 'CCC': ['P'], 'CCG': ['P'], 'CCU': ['P'], 'CGA': ['R'], 'CGC': ['R'], 'CGG': ['R'], 'CGU': ['R'], 'CUA': ['L'], 'CUC': ['L'], 'CUG': ['L'], 'CUU': ['L'], 'GAA': ['E'], 'GAC': ['D'], 'GAG': ['E'], 'GAU': ['D'], 'GCA': ['A'], 'GCC': ['A'], 'GCG': ['A'], 'GCU': ['A'], 'GGA': ['G'], 'GGC': ['G'], 'GGG': ['G'], 'GGU': ['G'], 'GUA': ['V'], 'GUC': ['V'], 'GUG': ['V'], 'GUU': ['V'], 'UAA': ['0'], 'UAC': ['Y'], 'UAG': ['0'], 'UAU': ['Y'], 'UCA': ['S'], 'UCC': ['S'], 'UCG': ['S'], 'UCU': ['S'], 'UGA': ['0'], 'UGC': ['C'], 'UGG': ['W'], 'UGU': ['C'], 'UUA': ['L'], 'UUC': ['F'], 'UUG': ['L'], 'UUU': ['F']}

mass_code = {'G': 57, 'A': 71, 'S': 87, 'P': 97, 'V': 99, 'T': 101, 'C': 103, 'I': 113, 'L': 113, 'N': 114, 'D': 115, 'K': 128, 'Q': 128, 'E': 129, 'M': 131, 'H': 137, 'F': 147, 'R': 156, 'Y': 163, 'W': 186}

spec_masses = {57: 'G', 71: 'A', 87: 'S', 97: 'P', 99: 'V', 101: 'T', 103: 'C', 113: 'L', 114: 'N', 115: 'D', 128: 'Q', 129: 'E', 131: 'M', 137: 'H', 147: 'F', 156: 'R', 163: 'Y', 186: 'W'}


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
    return sorted(lin_spec)


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



def expand_candidates(pep):
    expanded_peps = []
    for mass in spec_masses:
        expanded_pep = str(pep) + '-' + str(mass)
        expanded_peps.append(expanded_pep)
    return expanded_peps



def cand_seq(pep):
    seq = ''
    for mass in pep.split('-'):
        seq += spec_masses[int(mass)]
    return seq



def cyclopep_seq(spec):
    cand_peps = [str(mass) for mass in spec_masses]
    final_peps = []
    final_seq = []
    spec_int = [int(x) for x in spec]
    while len(cand_peps) > 0:
        for pep in cand_peps:
            total_mass = sum(int(x) for x in pep.split('-'))
            if total_mass not in spec_int:
                cand_peps.remove(pep)
            if total_mass in spec_int:
                if total_mass != max(spec_int):
                    check = all(mass in spec for mass in linear_spec(cand_seq(pep)))
                    if check:
                        for x in expand_candidates(pep):
                            cand_peps.append(x)
                        cand_peps.remove(pep)
                    else:
                        cand_peps.remove(pep)
                if total_mass == max(spec_int) and pep not in final_peps:
                    check = all(mass in spec for mass in linear_spec(cand_seq(pep)))
                    if check:
                        final_peps.append(pep)
                        final_seq.append(cand_seq(pep))
                        cand_peps.remove(pep)
                    else:
                        cand_peps.remove(pep)
            if total_mass > max(spec_int):
                pass
    answer = ''
    for pep in final_peps:
        answer += pep + ' '
    return answer


'''
Example
'''


print(cyclopep_seq([0, 113, 114, 128, 129, 227, 242, 242, 257, 355, 356, 370, 371, 484]))

# Output: '114-128-114-128 128-114-128-114 114-128-129-113 114-113-129-128
# 129-128-114-113 129-113-114-128 128-114-113-129 129-113-129-113 113-114-128-129
# 113-129-128-114 128-129-113-114 113-129-113-129'
