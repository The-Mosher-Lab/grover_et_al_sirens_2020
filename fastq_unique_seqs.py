#!/usr/bin/env python3

# Author: Jeffrey Grover
# Purpose: Output the unique sequences found in a .fastq file
# Created: 2020-03-04

from argparse import ArgumentParser
import gzip


def magic_open(input_file):
    if input_file.endswith('gz'):
        return gzip.open(input_file, 'rt')
    else:
        return open(input_file, 'r')


def fastq_count_seqs(input_fastq, min_len, max_len):
    profile_dict = {}
    with magic_open(input_fastq) as input_handle:
        n = 0
        for line in input_handle:
            n += 1
            if n == 2:
                seq = line.strip()
                if min_len <= len(seq) <= max_len:
                    if seq not in profile_dict:
                        profile_dict.update({seq: 0})
                    profile_dict[seq] += 1
            elif n == 4:
                n = 0
    return profile_dict


def output_profile(profile_dict):
    print('sequence', 'count', sep='\t')
    for sequence, count in profile_dict.items():
        print(sequence, count, sep='\t')


# Command line parser

def get_args():
    parser = ArgumentParser(
        description='Count unique sequences in a .fastq file.')
    parser.add_argument('fastq',
                        help='Input .fastq(.gz)',
                        metavar='FILE.fastq')
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
    profile = fastq_count_seqs(args.fastq, args.min_length, args.max_length)
    output_profile(profile)


if __name__ == '__main__':
    main(get_args())
