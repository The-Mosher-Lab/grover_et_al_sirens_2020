#!/usr/bin/env python3

# Author: Jeffrey Grover
# Purpose: Profile read legnths in a .bam alignment file
# Created: 2019-08-16

import pysam
from argparse import ArgumentParser
from os.path import exists
from sys import exit


def bam_length_profile(input_bam, min_len, max_len):
    profile_dict = {}
    for length in range(min_len, max_len + 1):
        profile_dict.update({length: 0})
    with pysam.AlignmentFile(input_bam, 'rb') as align_handle:
        for aln in align_handle.fetch():  # .fetch returns only mapped reads
            if min_len <= aln.query_length <= max_len:
                profile_dict[aln.query_length] += 1
    return profile_dict


def output_fastq_lengths(profile_dict):
    print('length', 'count', sep='\t')
    for length, count in profile_dict.items():
        print(length, count, sep='\t')


# Command line parser

def get_args():
    parser = ArgumentParser(
        description='Profile the mapped reads between a minimum and malixmum '
        'length in a .bam file.')
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

    profile = bam_length_profile(args.alignment, args.min_length, args.max_length)
    output_fastq_lengths(profile)


if __name__ == '__main__':
    main(get_args())
