#!/usr/bin/env python
# coding: utf-8

import sys
import random
from optparse import OptionParser

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Alphabet import IUPAC
from Bio.Geo.Record import out_block


def get_options():
    parser = OptionParser()
    parser.add_option("-i", "--input_filename", dest="input_filename",
                  help="input fasta filename", metavar="FILE")
    
    parser.add_option("-o", "--output_filename", dest="output_filename",
                  help="Output fasta filename", metavar="FILE")
    
    parser.add_option("-p", "--number_of_partitions", dest="partitions",
                  help="sampling percent (%)", type="int", default=0)
    
    parser.add_option("-s", "--number_of_sequences", dest="sequences",
                  help="number of subsequences ", type="int", default=0)
        
    return parser.parse_args()


def GetNumberOfSequences(input_handle, file_format):
    return len([None for record in SeqIO.parse(input_handle, file_format)])


def main():
    (options, args) = get_options()
    input_filename = options.input_filename
    if options.partitions == 0 and options.sequences == 0:
        sys.exit("error: you must set -p or -s.")
    file_format = "fasta"
    base_number_sequences = options.sequences
    input_handle = open(input_filename)
    number_remain_sequences = 0
    if options.partitions != 0 :
        number_sequences_in_file = GetNumberOfSequences(input_handle, file_format)
        base_number_sequences = number_sequences_in_file/options.partitions
        number_remain_sequences = number_sequences_in_file - (base_number_sequences*options.partitions)
        input_handle.close()
        input_handle = open(input_filename)
    
    
    file_id = 0
    sequence_count = 0
    output_handle = None
    number_sequences = base_number_sequences
    for record in SeqIO.parse(input_handle, file_format):
        if not output_handle :
            output_handle = open(options.output_filename + "_" + str(file_id), "w")
            number_sequences = base_number_sequences
            if file_id <  number_remain_sequences:
                number_sequences += 1
    
            file_id += 1
        if sequence_count >= number_sequences :
            sequence_count = 0
            output_handle.close()
            output_handle = open(options.output_filename + "_" + str(file_id), "w")
            number_sequences = base_number_sequences
            if file_id <  number_remain_sequences:
                number_sequences += 1
    
            file_id += 1
        SeqIO.write(record, output_handle, file_format)
        sequence_count += 1
        
    if output_handle : 
        output_handle.close()
    

if __name__ == '__main__':
    main()
        
