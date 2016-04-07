#!/usr/bin/env python
# file extract_primers.py

''' Extract Primers from sequence reads, for .fastq and .fasta files.
'''
__author__ = "Billy Fournier"
__license__ = "GPL"
__email__ = "billyfournier2000@yahoo.com"

from string import upper

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

header,mapping_file = process_mapping_file("mappingfile.txt")

print mapping_file[2]

if "LinkerPrimerSequence" in header:
    primer_ix = header.index("LinkerPrimerSequence")
print primer_ix

for primer in line[primer_ix].split(','):
    upper(primer).strip() 
