#!/usr/bin/env python

# var = {'A': 'A', 'N': 'ATCG'}
# lst = ['AN']

# var = {'A': 'A', 'N': 'ATCG', 'K': 'GT'}
# lst = ['ANK']

# var = {'A': 'A', 'W': 'AT', 'K': 'GT'}
# lst = ['AWK']

var = {'A': 'A', 'T': 'T', 'G': 'G', 'C': 'C', 'R': 'AG', 'Y': 'CT',
         'S': 'GC', 'W': 'AT', 'K': 'GT', 'M': 'AC', 'B': 'CGT',
         'D': 'AGT', 'H': 'ACT', 'V': 'ACG', 'N': 'ACGT'}

# var = {'A': 'A', 'V': 'ACG', 'K': 'GT'}
# lst = ['AKV']

# lst = ['CAGCMGCCGCGGTAA']
lst = ['TACNVGGGTATCTAATCC']


def get_string_count(sequence,nucleotide_variable):
    exp_str_total = 1
    for word in sequence:
        for ch in word:
            exp_str_total = exp_str_total * len(nucleotide_variable[ch])
    return exp_str_total

def create_perms(base_value,sequence_list):
    # Total sequence count needs to be divided by length
    # of base variations
    num_base = len(base_value)
    len_sl = len(sequence_list)
    cut_num = len_sl / num_base
    for index,items in enumerate(sequence_list):
        # index/cut_num
        base = index/cut_num
        items.append(var[ch][base])
    # print sequence_list, '\n'


sequence_list = []
for word in lst:
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
                    # print sequence_list
                create_perms(var[ch],sequence_list)
                print "\n***First Time Primer Finished***\n"


            # Adding "n" new sequences to sequence primer list
            else:
                print 'beginning of else:'
                to_add = len(var[ch]) * len(sequence_list) - len(sequence_list)
                for i in range( to_add ):
                    sequence_list.append(list( sequence_list[i] ))
                    # print sequence_list
                create_perms(var[ch],sequence_list)






for index,seq in enumerate(sequence_list):
    print index,seq
