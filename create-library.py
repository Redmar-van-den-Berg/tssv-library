#!/usr/bin/env python3

import argparse
import time
import os
import logging

from Bio import SeqIO

def create_library(name, chrom, begin, flank_size=20, end=None):
    """ Create a library line 

    chrom should be a Biopython Seq object of the relevant chromosome

    """
    if not end:
        end = begin + 1

    left_marker = chrom[begin - flank_size:begin].seq
    right_marker = chrom[end : end+flank_size].seq
    expected_allele = chrom[begin:end].seq

    return f'{left_marker}\t{right_marker}\t{expected_allele} 1 1'

def main(args):
    logging.debug(args)

    # Read the reference fasta (requires biopython >=1.52
    record_dict = SeqIO.index(args.reference, 'fasta')
    # Get the correct chromosome
    chrom = record_dict[args.chromosome]

    lib = create_library('name', chrom, args.begin, flank_size=args.flank_size)
    print(lib)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--reference', required=True,
                        help='Reference fasta file')
    parser.add_argument('--begin', required=True, type=int,
                        help='Beginning of region of interest (0 based)')
    parser.add_argument('--end', required=False, type=int,
                        help='End of region of interest (0 based)')
    parser.add_argument('--chromosome', default='chr7',
                        help='Chromosome for the region of interest')
    parser.add_argument('--flank-size', default=20, type=int,
                        help='Size of the flanking regions')

    arguments = parser.parse_args()
    # If no end is specified, 
    if not arguments.end:
        arguments.end = arguments.begin + 1
    main(arguments)
