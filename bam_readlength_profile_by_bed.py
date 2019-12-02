#!/usr/bin/env python3

# Author: Jeffrey Grover
# Purpose: Profile the read lengths which map to regions in a bed file
# Created: 2019-08-14
# Depends: pysam, samtools, python >= 3.6

import pysam
import gzip
from os.path import exists
from argparse import ArgumentParser
from sys import exit


def magic_open(input_file):
    if input_file.endswith('gz'):
        return gzip.open(input_file, 'rt')
    else:
        return open(input_file, 'r')


def bed_iter(input_file):
    with magic_open(input_file) as input_handle:
        for line in input_handle:
            entry = line.strip().split()
            yield {'chrom': entry[0],
                   'start': int(entry[1]),
                   'end': int(entry[2]),
                   'name': entry[3]}


def profile_reads_by_region(align_file, bed_iter, min_len, max_len):
    print('feature', *range(min_len, max_len + 1), sep='\t')
    for entry in bed_iter:
        profile = {}
        for length in range(min_len, max_len + 1):
            profile.update({length: 0})
        with pysam.AlignmentFile(align_file, 'rb') as align_handle:
            for aln in align_handle.fetch(entry['chrom'], entry['start'], entry['end']):
                read_length = aln.query_length
                if min_len <= read_length <= max_len:
                    profile[read_length] += 1
        print(*[entry['name']] + [profile[key] for key in profile], sep='\t')


# Command line parser

def get_args():
    parser = ArgumentParser(
        description='Profile the reads between a minimum and maximum length '
        'from regions defined by a .bed file.')
    parser.add_argument('alignment',
                        help='Input alignment file',
                        metavar='FILE.bam')
    parser.add_argument('-b', '--bed',
                        help='Input bed file containing regions of interest',
                        metavar='FILE.bed(.gz)')
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

    # Process files
    profile_reads_by_region(args.alignment, bed_iter(args.bed), args.min_length,
                            args.max_length)


if __name__ == '__main__':
    main(get_args())
