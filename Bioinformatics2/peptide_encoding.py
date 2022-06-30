'''
Peptide Encoding Solution
'''

'''
Given a DNA strand as a string, this function transcribes the strand its reverse complement
to determine what substrings in the DNA strand encode a specified peptide sequence. The function
also keeps track of the location of the peptide sequences within the entire genome.
'''



genetic_code = {'AAA': ['K'], 'AAC': ['N'], 'AAG': ['K'], 'AAU': ['N'], 'ACA': ['T'], 'ACC': ['T'], 'ACG': ['T'], 'ACU': ['T'], 'AGA': ['R'], 'AGC': ['S'], 'AGG': ['R'], 'AGU': ['S'], 'AUA': ['I'], 'AUC': ['I'], 'AUG': ['M'], 'AUU': ['I'], 'CAA': ['Q'], 'CAC': ['H'], 'CAG': ['Q'], 'CAU': ['H'], 'CCA': ['P'], 'CCC': ['P'], 'CCG': ['P'], 'CCU': ['P'], 'CGA': ['R'], 'CGC': ['R'], 'CGG': ['R'], 'CGU': ['R'], 'CUA': ['L'], 'CUC': ['L'], 'CUG': ['L'], 'CUU': ['L'], 'GAA': ['E'], 'GAC': ['D'], 'GAG': ['E'], 'GAU': ['D'], 'GCA': ['A'], 'GCC': ['A'], 'GCG': ['A'], 'GCU': ['A'], 'GGA': ['G'], 'GGC': ['G'], 'GGG': ['G'], 'GGU': ['G'], 'GUA': ['V'], 'GUC': ['V'], 'GUG': ['V'], 'GUU': ['V'], 'UAA': ['0'], 'UAC': ['Y'], 'UAG': ['0'], 'UAU': ['Y'], 'UCA': ['S'], 'UCC': ['S'], 'UCG': ['S'], 'UCU': ['S'], 'UGA': ['0'], 'UGC': ['C'], 'UGG': ['W'], 'UGU': ['C'], 'UUA': ['L'], 'UUC': ['F'], 'UUG': ['L'], 'UUU': ['F']}




def protein_trans(strand):
    code = genetic_code
    l = len(strand)
    protein = ''
    for i in range(l//3):
        codon = strand[i*3:(i*3)+3]
        if code[codon][0] != '0':
            protein += str(code[codon][0])
        else:
            break
    return protein



def reversecomplement(x):
    complement = ''
    for i in x:
        if i == 'A':
            complement += 'T'
        elif i == 'T':
            complement += 'A'
        elif i == 'C':
            complement += 'G'
        elif i == 'G':
            complement += 'C'
    return complement[::-1]



def transcribe(dna):
    rna = dna.replace('T', 'U')
    return rna


def peptide_code(dna, peptide):
    substrings = []
    substring_start_pos = 'Nucleotide start positions: '
    comp_dna = reversecomplement(dna)
    rna = transcribe(dna)
    comp_rna = transcribe(comp_dna)
    k = len(peptide)*3
    l = len(dna)
    for i in range(l - k + 1):
        string = rna[i:i+k]
        comp_string = comp_rna[i:i+k]
        if protein_trans(string) == peptide or protein_trans(comp_string) == peptide:
            substrings.append(dna[i:i+k])
        else:
            pass
    for i in range(l):
        if dna[i:i+k] in substrings:
            substring_start_pos += str(i) + ' '
    answer = ''
    for string in substrings:
        answer += string + ' '
    print(substring_start_pos)
    return answer


'''
Example
'''


print(peptide_code('GAAACTACCGAGATGCTACAATCGCTACGCGATCTG', 'AT'))

# Output: 'Nucleotide start positions: 14 23'
# 'GCTACA GCTACG'