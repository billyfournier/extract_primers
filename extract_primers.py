#!/usr/bin/env python

# USAGE:
# extract_primers.py <read1> <read2> <out1> <out2> <mappingfile>


from sys import argv
from string import upper
from re import compile

from skbio.parse.sequences import parse_fastq
from skbio.sequence import DNA
from skbio.format.sequences import format_fastq_record

def get_primers(mapping_file):
    """ Generates regular expression objects for forward and reverse primer
        from the mapping file.
    """
    iupac = {'A': 'A', 'T': 'T', 'G': 'G', 'C': 'C', 'R': '[AG]', 'Y': '[CT]',
             'S': '[GC]', 'W': '[AT]', 'K': '[GT]', 'M': '[AC]', 'B': '[CGT]',
             'D': '[AGT]', 'H': '[ACT]', 'V': '[ACG]', 'N': '[ACGT]'}
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
            # raw_forward_rc_primers.update([str(DNA(primer).rc()) for
            #                            primer in raw_forward_primers])
            raw_reverse_primers.update([upper(primer).strip() for
                                    primer in line[rev_primer_ix].split(',')])
            # raw_reverse_rc_primers.update([str(DNA(primer).rc()) for
            #                            primer in raw_reverse_primers])


    if not raw_forward_primers:
        raise ValueError(("No forward primers detected in mapping file."))
    if not raw_reverse_primers:
        raise ValueError(("No reverse primers detected in mapping file."))


    forward_primers = []
    reverse_primers = []

    for curr_primer in raw_forward_primers:
        forward_primers.append(''.join([iupac[symbol] for
                                                symbol in curr_primer]))
    for curr_primer in raw_reverse_primers:
        reverse_primers.append(''.join([iupac[symbol] for
                                                symbol in curr_primer]))

    forward_primers_re = compile("|".join([primer for primer in forward_primers]))
    reverse_primers_re = compile("|".join([primer for primer in reverse_primers]))

    return forward_primers_re, reverse_primers_re

def remove_primers(input_fastq, output_fastq,primers):
    count = 0
    # USING regex list (Time 11m4)
    with open(input_fastq) as read, open(output_fastq, "w") as out_seqs:
        for label,seq,qual in parse_fastq(read):
            start_slice = 0
            if primers.search(seq):
                start_slice = int(primers.search(seq).end())
            curr_seq = seq[start_slice:]
            curr_qual = qual[start_slice:]
            if start_slice > 0:
                formatted_fastq_line = format_fastq_record(label, curr_seq, curr_qual)
                # print ("%s" % (formatted_fastq_line))
                out_seqs.write("%s" % (formatted_fastq_line))

################################################################################
####                                                                        ####
################################################################################
f_fastq1 = argv[1]
r_fastq2 = argv[2]
out_seqs1 = argv[3]
out_seqs2 = argv[4]
map_file = argv[5]

forward_primers, reverse_primers = get_primers(map_file)
remove_primers(f_fastq1,out_seqs1,forward_primers)
remove_primers(r_fastq2,out_seqs2,reverse_primers)
