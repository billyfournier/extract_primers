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

from skbio.parse.sequences import parse_fastq
from skbio.sequence import DNA
from skbio.format.sequences import format_fastq_record

from qiime.check_id_map import process_id_map
from qiime.split_libraries_fastq import (check_header_match_pre180,
                                         check_header_match_180_or_later)
from qiime.parse import is_casava_v180_or_later
from qiime.pycogent_backports.fastq import FastqParseError

def get_mappingfile_header(mapping_data):
    fin = open(mapping_data)
    header = fin.readline()
    return header

# unit test. Be sure mappingfile
#print get_mappingfile_header("mappingfile.txt")
def get_mapping_data(mapping_file):
    header, mapping_data, run_description, errors, warnings =\
        process_id_map(mapping_file)
    return mapping_data



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
        raw_forward_rc_primers.update([str(DNA(primer).rc()) for
                                       primer in raw_forward_primers])
        raw_reverse_primers.update([upper(primer).strip() for
                                    primer in line[rev_primer_ix].split(',')])
        raw_reverse_rc_primers.update([str(DNA(primer).rc()) for
                                       primer in raw_reverse_primers])

    if not raw_forward_primers:
        raise ValueError(("No forward primers detected in mapping file."))
    if not raw_reverse_primers:
        raise ValueError(("No reverse primers detected in mapping file."))

    # Finding the forward primers, or rc of reverse primers indicates forward
    # read. Finding the reverse primer, or rc of the forward primers, indicates
    # the reverse read, so these sets are merged.
    raw_forward_primers.update(raw_reverse_rc_primers)
    raw_reverse_primers.update(raw_forward_rc_primers)

    forward_primers = []
    reverse_primers = []
    for curr_primer in raw_forward_primers:
        forward_primers.append(compile(''.join([iupac[symbol] for
                                                symbol in curr_primer])))
    for curr_primer in raw_reverse_primers:
        reverse_primers.append(compile(''.join([iupac[symbol] for
                                                symbol in curr_primer])))

    return forward_primers, reverse_primers



header = get_mappingfile_header("mappingfile.txt")
mapping_data = get_mapping_data("mappingfile.txt")

print get_primers(header,mapping_data)









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
