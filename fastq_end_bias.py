#!/usr/bin/env python3

# Author: Jeffrey Grover
# Purpose: Profile nucleotide bias at ends of reads from a .fastq file
# Created: 2019-08-09

import gzip
from argparse import ArgumentParser


def magic_open(input_file):
    if input_file.endswith('gz'):
        return gzip.open(input_file, 'rt')
    else:
        return open(input_file, 'r')


def fastq_yield_seqs(input_fastq):
    with magic_open(input_fastq) as input_handle:
        i = 0
        for line in input_handle:
            i += 1
            if i == 2:
                yield line.strip()
            elif i == 4:
                i = 0


def end_bias(fastq_seqs, min_length, max_length):
    bias_dict = {'5_prime': {'A': 0, 'T': 0, 'C': 0, 'G': 0},
                 '3_prime': {'A': 0, 'T': 0, 'C': 0, 'G': 0}}
    for sequence in fastq_seqs:
        if min_length <= len(sequence) <= max_length:
            bias_dict['5_prime'][sequence[0]] += 1
            bias_dict['3_prime'][sequence[-1]] += 1
    return bias_dict


# Command line parser

def get_args():
    parser = ArgumentParser(
        description='Profile 3\' and 5\' nucleotide bias from .fastq file')
    parser.add_argument('fastq',
                        help='Input .fastq, may be gzipped',
                        metavar='FILE.fastq(.gz)')
    parser.add_argument('-n', '--min_length',
                        help='Minimum length of reads to profile (default=0)',
                        default=0,
                        type=int,
                        metavar='INT')
    parser.add_argument('-m', '--max_length',
                        help='Maximum length of reads to profile (default=150)',
                        default=150,
                        type=int,
                        metavar='INT')
    return parser.parse_args()


# Main function entry point

def main(args):
    bias_dict = end_bias(
        fastq_yield_seqs(args.fastq), args.min_length, args.max_length
    )
    print('end', 'A', 'T', 'C', 'G', sep=',')
    for end, freq_dict in bias_dict.items():
        line = [end]
        for base, count in freq_dict.items():
            line.append(str(count))
        print(','.join(line))


if __name__ == '__main__':
    main(get_args())
