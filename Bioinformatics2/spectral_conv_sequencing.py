'''
Spectral Convolution Branch and Bound Algorithm Solution
'''

'''
Many antibodies are cyclic peptides that do not conform to the standards present in the central
dogma of molecular biology. This includes utilizing amino acids outside of the standard set of
20 that are commonly discussed in translation. To account for this, the following algorithm
forms a convolution spectrum and utilizes the most commonly occurring masses in order to begin
building candidate peptides to identify a peptide sequence from a given experimental spectrum.
Additionally, the list of candidate peptides is trimmed after each branching step so only the
best performing candidates are kept. Due to memory and time constraints, a break point was inserted
once a certain number of worse performing candidates were compared to the best performing peptide
sequence.
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
            return mass_list[:i]
    return mass_list



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



def conv_pep_seq(spectrum, m, n):
    mass_list = spectral_conv(spectrum, m)
    leaderboard = [str(mass) for mass in mass_list]
    lead_pep = ''
    best_peps = []
    x = 0
    spec_int = [int(x) for x in spectrum]
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


print(conv_pep_seq(['0', '113', '114', '128', '129', '242', '242', '257', '370', '371', '484'], 10, 10))


# Output: 'The best scoring sequence is: 113-129-128-114
# Alternate peptide solutions: ['113-129-128-114', '114-128-129-113']'
