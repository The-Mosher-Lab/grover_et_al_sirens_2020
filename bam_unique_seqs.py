#!/usr/bin/env python3

# Author: Jeffrey Grover
# Purpose: Output the unique aligned sequences found in a .bam file
# Created: 2020-03-04

import pysam
from argparse import ArgumentParser
from os.path import exists
from sys import exit


def bam_count_seqs(input_bam, min_len, max_len):
    profile_dict = {}
    with pysam.AlignmentFile(input_bam, 'rb') as align_handle:
        for aln in align_handle.fetch():
            if min_len <= aln.query_length <= max_len:
                if aln.query_sequence not in profile_dict:
                    profile_dict.update({aln.query_sequence: 0})
                profile_dict[aln.query_sequence] += 1
    return profile_dict


def output_aligned_profile(profile_dict):
    print('sequence', 'count', sep='\t')
    for sequence, count in profile_dict.items():
        print(sequence, count, sep='\t')


# Command line parser

def get_args():
    parser = ArgumentParser(
        description='Count unique sequences that have been aligned to a reference in a .bam file')
    parser.add_argument('alignment',
                        help='Input alignment file',
                        metavar='FILE.bam')
    parser.add_argument('-n', '--min_length',
                        help='Minimum length of reads to profile',
                        type=int,
                        metavar='INT')
    parser.add_argument('-m', '--max_length',
                        help='Maximum length of reads to profile',
                        type=int,
                        metavar='INT')
    return parser.parse_args()


# Main function entry point

def main(args):

    # Check that alignment file is a bam
    if not args.alignment.endswith('.bam'):
        exit('Error: Alignment must be in .bam format to enable random access.')

    # Create .bai index if needed
    if not exists(args.alignment + '.bai'):
        pysam.index(args.alignment)

    profile = bam_count_seqs(args.alignment, args.min_length, args.max_length)
    output_aligned_profile(profile)


if __name__ == '__main__':
    main(get_args())
