#!/usr/bin/env python
# file extract_primers.py

''' Extract Primers from sequence reads, for .fastq and .fasta files.
'''
__author__ = "Billy Fournier"
__license__ = "GPL"
__email__ = "billyfournier2000@yahoo.com"

from string import upper

def process_mapping_file(mapping_file):
    with open(mapping_file) as mappingfile:
        header = mappingfile.readline().split()

        mapping_data = []
        mappingfile.readline() ### Could cause issues. ###
        for line in mappingfile:
            mapping_data.append(list(line.split()))
    return header, mapping_data

header,mapping_data = process_mapping_file("mappingfile.txt")

print mapping_data[2]

if "LinkerPrimerSequence" in header:
    primer_ix = header.index("LinkerPrimerSequence")
print primer_ix

for primer in line[primer_ix].split(','):
    upper(primer).strip()
