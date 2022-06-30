'''
Leaderboard Branch and Bound Algorithm Solution
'''

'''
This algorithm utilizes a branch and bound solution to build candidate peptides one amino acid at a time
from a given mass list and compare them for accuracy to a given experimental spectrum. At each step,
the 'leaderboard' of the possible peptide is trimmed to prevent the possible list of candidates from
growing too large while still accounting for errors in the experimental spectrum. The final peptide
returned is one which scores optimally with the given experimental spectrum. 
'''


mass_code = {'G': 57, 'A': 71, 'S': 87, 'P': 97, 'V': 99, 'T': 101, 'C': 103, 'I': 113, 'L': 113, 'N': 114, 'D': 115, 'K': 128, 'Q': 128, 'E': 129, 'M': 131, 'H': 137, 'F': 147, 'R': 156, 'Y': 163, 'W': 186}

spec_masses = {57: 'G', 71: 'A', 87: 'S', 97: 'P', 99: 'V', 101: 'T', 103: 'C', 113: 'L', 114: 'N', 115: 'D', 128: 'Q', 129: 'E', 131: 'M', 137: 'H', 147: 'F', 156: 'R', 163: 'Y', 186: 'W'}



def expand_leaderboard(leaderboard, mass_list):
    expanded_leaderboard = []
    for pep in leaderboard:
        for mass in mass_list:
            expanded_pep = str(pep) + '-' + str(mass)
            expanded_leaderboard.append(expanded_pep)
    return expanded_leaderboard



def conv_cyclo_spec(peptide):
    masses = [int(x) for x in peptide.split('-')]
    l = len(masses)
    pref_mass = [0]
    for i in range(l):
        pref_mass.append(pref_mass[i] + masses[i])
    cyclo_spec = [0]
    total_mass = pref_mass[l]
    for i in range(l):
        for j in range(i+1, l+1):
            cyclo_spec.append(pref_mass[j] - pref_mass[i])
            if i > 0 and j < l:
                cyclo_spec.append(total_mass - (pref_mass[j] - pref_mass[i]))
    return sorted(cyclo_spec)



def conv_lin_spec(peptide):
    masses = [int(x) for x in peptide.split('-')]
    l = len(masses)
    pref_mass = [0]
    for i in range(l):
        pref_mass.append(pref_mass[i] + masses[i])
    lin_spec = [0]
    for i in range(l):
        for j in range(i+1, l+1):
            lin_spec.append(pref_mass[j] - pref_mass[i])
    return sorted(lin_spec)



def linear_score(peptide, spectrum):
    if len(peptide) == 0:
        return 0
    theory_spec = conv_lin_spec(peptide)
    # print(theory_spec)
    score = 0
    x = 0
    for i in range(len(theory_spec)):
        for j in range(x, len(spectrum)):
            if str(spectrum[j]) == str(theory_spec[i]):
                score += 1
                x = j + 1
                break
    return score


def leaderboard_trim(leaderboard, spectrum, n):
    l = len(leaderboard)
    if l <= int(n):
        return leaderboard
    scores_dict = {}
    for i in range(l):
        peptide = leaderboard[i]
        scores_dict.setdefault(peptide, 0)
        scores_dict[peptide] = linear_score(peptide, spectrum)
    sorted_scores = dict(sorted(scores_dict.items(), key=lambda x: x[1], reverse=True))
    leaderboard = list(sorted_scores.keys())
    linear_scores = list(sorted_scores.values())
    k = len(linear_scores)
    cutoff = linear_scores[n-1]
    for i in range(n, k):
        if linear_scores[i] < cutoff:
            return leaderboard[:i]
    return leaderboard



def lead_cyclo_seq(spectrum, n):
    mass_list = spec_masses
    leaderboard = [str(mass) for mass in mass_list]
    lead_pep = ''
    best_peps = []
    spec_int = [int(x) for x in spectrum]
    x = 0
    while len(leaderboard) < 10000:
        leaderboard = expand_leaderboard(leaderboard, mass_list)
        for pep in leaderboard:
            total_mass = sum(int(x) for x in pep.split('-'))
            if total_mass > max(spec_int):
                leaderboard.remove(pep)
            elif total_mass == max(spec_int):
                if linear_score(pep, spectrum) > linear_score(lead_pep, spectrum):
                    best_peps = [pep]
                    lead_pep = pep
                elif linear_score(pep, spectrum) == linear_score(lead_pep, spectrum):
                    best_peps.append(pep)
                else:
                    x += 1
        leaderboard = leaderboard_trim(leaderboard, spectrum, n)
        if x > 100:
            break
    print('The best scoring sequence is: ' + lead_pep)
    return 'Alternate peptide solutions: ' + str(best_peps)


'''
Example
'''


print(lead_cyclo_seq(['0', '113', '114', '128', '129', '242', '242', '257', '370', '371', '484'], 100))


# Output: 'The best scoring sequence is: 113-129-128-114
# Alternate peptide solutions: ['113-129-128-114', '114-128-129-113', '113-129-128-57-57',
# '113-129-57-71-114', '113-129-71-57-114', '57-57-128-129-113', '114-57-71-129-113', '114-71-57-129-113']'
