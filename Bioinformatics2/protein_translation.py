'''
Protein Translation Solution
'''

'''
Given a RNA sequence as a string, this function returns the amino acid sequence that results
from the RNA strand using the genetic code to build a protein one codon at a time. This
function also accounts for stop codons by terminating translation once the sequence for
a stop codon is reached.
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
    return 'Amino Acid Sequence: ' + protein


'''
Example
'''


print(protein_trans('CCCAGGACUGAGAUCAAU'))

# Output: 'Amino Acid Sequence: PRTEIN'