#!/usr/bin/env python

# var = {'A': 'A', 'N': 'ATCG'}
# lst = ['AN']

var = {'A': 'A', 'N': 'ATCG', 'K': 'GT'}
lst = ['ANK']



def get_string_count(sequence,nucleotide_variable):
    exp_str_total = 1
    for word in sequence:
        for ch in word:
            exp_str_total = exp_str_total * len(nucleotide_variable[ch])
    return exp_str_total
# print get_string_count(lst,var)


sequence_list = []
sequence = []

for word in lst:
    sequence_list = [[]]

    for ch in word:
        # print "ch in word:"
        to_add = (len(var[ch]) * len(sequence_list)) - len(sequence_list)
        # print "to_add is: "
        # print to_add

# spot is not variable, so copy it to all sequence copies
        if len(var[ch]) is 1:
            # print len(var[ch])
            # print "IF"
            for string in sequence_list:
                string.append(ch)
            # print sequence_list
# spot is variable, so must create new sequence primers
        else:
            # print "ELSE"
            # print sequence_list
            # First time primers are added. **Edge case**
            if len(sequence_list) is 1:
                # print "if"
                # print sequence_list
                for i in range(1,4):
                    sequence_list.append([])
                    # print ''.join(sequence_list[0])
                    sequence_list[i].append(''.join(sequence_list[0]))
                    # print sequence_list

            # create n copies of sequence primer list
            else:
                for i in range( len(var[ch]) - 1 ):
                    for j in range(len(sequence_list)):
                        sequence_list.append(sequence_list[j])
                # print "else"
                # print sequence_list


            # print "\n\n"
            # print sequence_list
            # for index,items in enumerate(sequence_list):
            #     print index, items
            #     good_index = index % len(var[ch])
            #     print var[ch][good_index]
            #     print sequence_list
            #
            #     items.append(var[ch][good_index])
            #     print items
            #     print '\n'
            # print sequence_list

            # print "\n\n"
            # print "starting with", sequence_list
            # for index,item in enumerate(var[ch]):
            #     print "index is: ", index, "\titem is: ", item
            #     for index2, seq in enumerate(sequence_list):
            #         # print "len(item) is: ", len(item)
            #         # print "index2+1 is: ", index2+1
            #         good_index = index2 / len(var[ch])
            #         print "length: ",len(var[ch])
            #         # print "good_index is: " , good_index
            #         # print items[good_index]
            #         if good_index is index:
            #             seq.append(item)
            #             print sequence_list
            # print "ending with", sequence_list



# print sequence_list
