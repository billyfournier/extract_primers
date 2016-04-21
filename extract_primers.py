#!/usr/bin/env python
# file extract_primers.py
''' Extract Primers from sequence reads, for .fastq and .fasta files.
'''
__author__ = "Billy Fournier"
__license__ = "GPL"
__email__ = "billyfournier2000@yahoo.com"
__credits__ = "William Walters"

from string import upper
from itertools import product
from optparse import OptionParser
from skbio.sequence import DNA

parser = OptionParser()
parser.add_option("-f", "--forward_read", dest="forward_file",
                        help="(default) read1.fastq", default="read1.fastq")
parser.add_option("-r", "--reverse_read", dest="reverse_file",
                        help="(default) read2.fastq", default="read2.fastq")
parser.add_option("-o", "--forward_out", dest="forward_out",
                        help="name of forward output file", default="out1.fastq")
parser.add_option("-O", "--reverse_out", dest="reverse_out",
                        help="name of reverse output file", default="out2.fastq")
parser.add_option("-m", "--mapping_file", dest="mapping_file",
                        default="mappingfile.txt", help="mappingfile for fastq")

(options, args) = parser.parse_args()


# if options.input_file is None:
# 	print "ERROR: flag (-i) is required. "
# 	print "       This is the PATH to your read.fastq file"
# 	exit(-1)
if options.mapping_file is None:
	print "ERROR: flag (-m) is required. "
	print "       This is the PATH to your mappingfile.txt"
	exit(-1)
# if options.output_file is None:
# 	print "ERROR: flag (-o) is required. "
# 	print "       This is the name of the output"
# 	exit(-1)

""" Creates all combinations of a primer sequence with a known variable base/s.

EXAMPLE:
A primer "ATCGN" has 4 variations at the last position so 4 bases can be
represented. Those primers would be:
"ATCGA"     "ATCGT"     "ATCGC"     "ATCGG"
"""

def gen_primer_list(mapping_primer):
    var = {'A': 'A', 'T': 'T', 'G': 'G', 'C': 'C', 'R': 'AG', 'Y': 'CT',
             'S': 'GC', 'W': 'AT', 'K': 'GT', 'M': 'AC', 'B': 'CGT',
             'D': 'AGT', 'H': 'ACT', 'V': 'ACG', 'N': 'ACGT'}
    if type(mapping_primer) is list:
        poss = [var[c] for c in mapping_primer[0]]
    if type(mapping_primer) is str:
        poss = [var[c] for c in mapping_primer]
    return list(product(*poss))

"""Returns the reverse complement of a nucleotide sequence

The string sequence's complement is generated and then
the order of the compment is reversed.

Example:
sequence "AATTCC" is used to create its complement:
complement "TTAAGG", then this complement is reversed:
reverse complement "GGAATT".
"""
complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
def reverse_complement(seq):
    bases = list(seq)
    bases = reversed([complement.get(base,base) for base in bases])
    bases = list(bases)
    return bases


def get_primers(mapping_file):
    """ Generates and returns all forward and reverse primer combinations
    from the mapping file.
    mapping_data:  the sequenced data's mapping file
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
    raw_forward_rc_primers = set([])
    raw_reverse_primers = set([])
    raw_reverse_rc_primers = set([])

    for line in mapping_data:
        # Split on commas to handle pool of primers
        if len(line) >= rev_primer_ix:
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

    forward_primers = set([])
    reverse_primers = set([])
    for sequence in raw_forward_primers:
        forward_primers |= set(gen_primer_list(sequence))
    for sequence in raw_reverse_primers:
        reverse_primers |= set(gen_primer_list(sequence))

    return forward_primers, reverse_primers


# sequence identifier beginning format
# @<instrument>:<run number>:
def remove_primers(read_file,out_name,primers,tolerance):
    is_sequence = False
    seq_passed = False
    seq_identifier = ''
    seq_data = ''
    qual_data = ''
    # print index_file
    with open(read_file) as input:
        tmp_str = input.readline()
        index = tmp_str.find(':')
        instrument_id = tmp_str[0:index]

    with open(read_file) as input, open(out_name,"w") as output:
        for index, line in enumerate(input):
            # write qual data & Marker for sequence data

            if instrument_id in line:
                is_sequence = True
                seq_identifier = line
                if seq_passed:
                    qual_data = qual_data[new_index:]
                    output.write('+\n')

                    # print ''.join(qual_data.split())
                    output.write(qual_data)
                seq_passed = False

            # Marker for quality data and performs primer check
            elif line.split() == ['+']:
                is_sequence = False
                for primer in forward_primers:
                    new_index = seq_data.rfind(''.join(primer))
                    if new_index is not -1:
                        new_index = new_index + len(primer) - 1
                        seq_data = seq_data[new_index:]
                        if len(seq_data) > tolerance:
                            # print ''.join(seq_identifier.split())
                            # print ''.join(seq_data.split())
                            output.write(seq_identifier)
                            output.write(seq_data)
                            seq_passed = True
                            break

            # Populate seq_data with sequence string
            elif is_sequence:
                seq_data = seq_data + line
                qual_data = ''

            # Populates qual_data & clears seq_data
            elif not is_sequence:
                qual_data = qual_data + line
                seq_data = ''

        if seq_passed:
            qual_data = qual_data[new_index:]
            # print "+"
            # print ''.join(qual_data.split())
            output.write('+\n')
            output.write(qual_data)

mapping_file = options.mapping_file
forward_read = options.forward_file
reverse_read = options.reverse_file
forward_out = options.forward_out
reverse_out = options.reverse_out

forward_primers, reverse_primers = get_primers(mapping_file)
remove_primers(forward_read, forward_out, forward_primers, 50)
remove_primers(reverse_read, reverse_out, reverse_primers, 50)
