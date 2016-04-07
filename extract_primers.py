#!/usr/bin/env python
# file extract_primers.py
from __future__ import division

''' Extract Primers from sequence reads, for .fastq and .fasta files.
'''
__author__ = "Billy Fournier"
__license__ = "GPL"
__email__ = "billyfournier2000@yahoo.com"



import numpy as np

from string import upper
from itertools import izip, cycle
from os.path import join
from os import rename
from re import compile


def get_mappingfile_header(mapping_data):
    fin = open(mapping_data)
    header = fin.readline()
    return header


def process_mapping_file(mapping_data):
    with open(mapping_data) as mappingfile:
        header = mappingfile.readline()
        header = header.split()

        mapping_list = []
        mapping_element = []
        mappingfile.readline()
        for line in mappingfile:
            str = line
            str = str.split()
            for word in str:
                mapping_element.append(word)
            mapping_list.append(mapping_element)
            mapping_element = []
    return header, mapping_list


def get_primers(header,
                mapping_data):
    """ Returns lists of forward/reverse primer regular expression generators
    header:  list of strings of header data.
    mapping_data:  list of lists of mapping data
    Will raise error if either the LinkerPrimerSequence or ReversePrimer fields
        are not present
    """

    if "LinkerPrimerSequence" in header:
        primer_ix = header.index("LinkerPrimerSequence")
        print primer_ix
    else:
        raise IndexError(
            ("Mapping file is missing LinkerPrimerSequence field."))
    if "ReversePrimer" in header:
        rev_primer_ix = header.index("ReversePrimer")
    else:
        raise IndexError(("Mapping file is missing ReversePrimer field."))

    iupac = {'A': 'A', 'T': 'T', 'G': 'G', 'C': 'C', 'R': '[AG]', 'Y': '[CT]',
             'S': '[GC]', 'W': '[AT]', 'K': '[GT]', 'M': '[AC]', 'B': '[CGT]',
             'D': '[AGT]', 'H': '[ACT]', 'V': '[ACG]', 'N': '[ACGT]'}

    raw_forward_primers = set([])
    raw_forward_rc_primers = set([])
    raw_reverse_primers = set([])
    raw_reverse_rc_primers = set([])

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

    # Finding the forward primers, or rc of reverse primers indicates forward
    # read. Finding the reverse primer, or rc of the forward primers, indicates
    # the reverse read, so these sets are merged.
    # raw_forward_primers.update(raw_reverse_rc_primers)
    # raw_reverse_primers.update(raw_forward_rc_primers)

    print raw_forward_primers
    forward_primers = []
    reverse_primers = []
    
	for curr_primer in raw_forward_primers:
        forward_primers.append(compile(''.join([iupac[symbol] for
                                                symbol in curr_primer])))
    for curr_primer in raw_reverse_primers:
        reverse_primers.append(compile(''.join([iupac[symbol] for
                                                symbol in curr_primer])))
    print forward_primers

    return forward_primers, reverse_primers


header,mapping_data = process_mapping_file("mappingfile.txt")

forward_primers, reverse_primers = get_primers(header,mapping_data)
print forward_primers









#largestPrimer = max(len(str) for str in primers)

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
