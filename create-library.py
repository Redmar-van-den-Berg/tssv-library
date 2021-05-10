#!/usr/bin/env python3

import argparse

from cyvcf2 import VCF
from Bio import SeqIO

def create_library(name, chrom, begin, end=None, flank_size=20):
    """ Create a library line 

    chrom should be a Biopython Seq object of the relevant chromosome

    """
    if not end:
        end = begin + 1

    left_marker = chrom[begin - flank_size:begin].seq
    right_marker = chrom[end : end+flank_size].seq
    expected_allele = chrom[begin:end].seq

    return f'{name}\t{left_marker}\t{right_marker}\t{expected_allele} 1 1'

def main(args):
    # Read the reference fasta (requires biopython >=1.52
    record_dict = SeqIO.index(args.reference, 'fasta')

    # Open the VCF
    for variant in VCF(args.vcf):
        name = f'{variant.CHROM}:{variant.POS}'
        chrom = record_dict[variant.CHROM]
        # Get the size of the event
        size = max (len(x) for x in [variant.REF]+variant.ALT)
        # Skip if the event is too large
        if size > args.max_size:
            continue
        # VCF's are 1-based
        begin = variant.POS - 1
        end = begin + size
        lib = create_library(name, chrom, begin, end, args.flank_size)
        print(lib)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--reference', required=True,
                        help='Reference fasta file')
    parser.add_argument('--vcf', required=True,
                        help='Create a library for the variants in the VCF')
    parser.add_argument('--flank-size', default=20, type=int,
                        help='Size of the flanking regions')
    parser.add_argument('--max-size', default=20, type=int,
                        help='Maximum size of an indel to verify')

    arguments = parser.parse_args()
    main(arguments)
