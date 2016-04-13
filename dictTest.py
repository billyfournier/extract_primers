#!/usr/bin/env python
# lst = ['CAGCMGCCGCGGTAA']
lst = ['TACNVGGGTATCTAATCC']

#
#
#

# Method using itertools
from itertools import product
word = "AWK"
# word = ['AWK']

''' Produces a list of primers. 

    ARGUMENTS:
        mapping_primer -- This is the raw primer from the mapping file.
    RETURN:
        Return type is list
'''
def gen_primer_list(mapping_primer):
    var = {'A': 'A', 'T': 'T', 'G': 'G', 'C': 'C', 'R': 'AG', 'Y': 'CT',
             'S': 'GC', 'W': 'AT', 'K': 'GT', 'M': 'AC', 'B': 'CGT',
             'D': 'AGT', 'H': 'ACT', 'V': 'ACG', 'N': 'ACGT'}
    if type(mapping_primer) is list:
        poss = [var[c] for c in mapping_primer[0]]
    if type(mapping_primer) is str:
        poss = [var[c] for c in mapping_primer]
    return list(product(*poss))

print gen_primer_list(word)
# >>> [''.join(p) for p in product(*poss)]
# ['AAG', 'AAT', 'ATG', 'ATT']


#
#
#

# Long hand method
def create_perms(base_value,sequence_list):
    # Total sequence count needs to be divided by length
    # of base variations
    num_base = len(base_value)
    len_sl = len(sequence_list)
    cut_num = len_sl / num_base
    for index,items in enumerate(sequence_list):
        base = index/cut_num
        items.append(base_value[base])

def get_primer_list(primer):
    var = {'A': 'A', 'T': 'T', 'G': 'G', 'C': 'C', 'R': 'AG', 'Y': 'CT',
             'S': 'GC', 'W': 'AT', 'K': 'GT', 'M': 'AC', 'B': 'CGT',
             'D': 'AGT', 'H': 'ACT', 'V': 'ACG', 'N': 'ACGT'}
    sequence_list = []
    for word in primer:
        sequence_list = [[]]

        for ch in word:
            # Base is not variable so append it to all seq
            if len(var[ch]) is 1:
                for seq in sequence_list:
                    seq.append(ch)

            else:    # Base is variable so...
                # First time primers are added. **Edge case**
                if len(sequence_list) is 1:
                    to_add = len(var[ch]) - 1
                    for i in range(to_add):
                        sequence_list.append(list( sequence_list[0] ))
                    create_perms(var[ch],sequence_list)
                    print "\n***First Time Primer Finished***\n"


                # Adding "n" new sequences to sequence primer list
                else:
                    print 'beginning of else:'
                    to_add = len(var[ch]) * len(sequence_list) - len(sequence_list)
                    for i in range( to_add ):
                        sequence_list.append(list( sequence_list[i] ))
                    create_perms(var[ch],sequence_list)
    return sequence_list

# >>> from itertools import product
# >>> var = {'A': 'A', 'W': 'AT', 'K': 'GT'}
# >>> word = "AWK"
# >>> poss = [var[c] for c in word]
# >>> poss
# ['A', 'AT', 'GT']
# >>> list(product(*poss))
# [('A', 'A', 'G'), ('A', 'A', 'T'), ('A', 'T', 'G'), ('A', 'T', 'T')]

# >>> [''.join(p) for p in product(*poss)]
# ['AAG', 'AAT', 'ATG', 'ATT']

#
# seqs = get_primer_list(lst)
# for index,seq in enumerate(seqs):
#     print index,seq, "\n"
