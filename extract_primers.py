#!/usr/bin/env python
# file extract_primers.py
''' Extract Primers from sequence reads, for .fastq and .fasta files.
'''
__author__ = "Billy Fournier"
__license__ = "GPL"
__email__ = "billyfournier2000@yahoo.com"

from string import upper
from itertools import product


def gen_primer_list(mapping_primer):
    var = {'A': 'A', 'T': 'T', 'G': 'G', 'C': 'C', 'R': 'AG', 'Y': 'CT',
             'S': 'GC', 'W': 'AT', 'K': 'GT', 'M': 'AC', 'B': 'CGT',
             'D': 'AGT', 'H': 'ACT', 'V': 'ACG', 'N': 'ACGT'}
    if type(mapping_primer) is list:
        poss = [var[c] for c in mapping_primer[0]]
    if type(mapping_primer) is str:
        poss = [var[c] for c in mapping_primer]
    return list(product(*poss))


def get_primers(mapping_file):
    """ Returns lists of forward/reverse primer regular expression generators
    header:  list of strings of header data.
    mapping_data:  list of lists of mapping data
    Will raise error if either the LinkerPrimerSequence or ReversePrimer fields
        are not present
    """
    # processing mapping_file into header and mapping_data
    with open(mapping_file) as mappingfile:
        header = mappingfile.readline().split()
        mapping_data = []
        mappingfile.readline() ### Could cause issues. ###
        for line in mappingfile:
            mapping_data.append(list(line.split()))


    if "LinkerPrimerSequence" in header:
        primer_ix = header.index("LinkerPrimerSequence") # ix = index
    else:
        raise IndexError(
            ("Mapping file is missing LinkerPrimerSequence field."))
    if "ReversePrimer" in header:
        rev_primer_ix = header.index("ReversePrimer")
    else:
        raise IndexError(("Mapping file is missing ReversePrimer field."))

    raw_forward_primers = set([])
    raw_reverse_primers = set([])
    for line in mapping_data:
        # Split on commas to handle pool of primers
        raw_forward_primers.update([upper(primer).strip() for
                                    primer in line[primer_ix].split(',')])
        raw_reverse_primers.update([upper(primer).strip() for
                                    primer in line[rev_primer_ix].split(',')])

    if not raw_forward_primers:
        raise ValueError(("No forward primers detected in mapping file."))
    if not raw_reverse_primers:
        raise ValueError(("No reverse primers detected in mapping file."))

    forward_primers = set([])
    reverse_primers = set([])
    for sequence in raw_forward_primers:
        forward_primers |= set(gen_primer_list(sequence))
    for sequence in raw_reverse_primers:
        reverse_primers |= set(gen_primer_list(sequence))

    return forward_primers, reverse_primers



forward_primers, reverse_primers = get_primers("mappingfile.txt")


for seq in forward_primers:
    print seq
# for seq in reverse_primers:
#     print seq

def removePrimers(inFile,outFile,primerList,tolerance):
    items = ["@","+"]
    with open(inFile) as input, open(outFile,"w") as output:
        for line in input:
            if line[0] in items:
                output.write(line)
            else:
                str = line
                indexList = []
                for primer in primers:
                    length = len(primer)
                    for i in range(len(str)-length):
                        window = str[i:i+length]
                        if hammingDistance(window,primer) <= tolerance:
                            indexList.append(i+length)
                if len(indexList) is 2:
                    output.write(str[indexList[0]:indexList[1]] + "\n")

### unit test
'''
primers = ["gattag","taggat"]
removePrimers("test","testout",primers,0)
'''
