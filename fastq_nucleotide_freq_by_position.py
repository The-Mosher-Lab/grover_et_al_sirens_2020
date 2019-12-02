#!/usr/bin/env python3

# Author: Jeffrey Grover
# Purpose: Profile nucleotide bias by position from reads in a .fastq file
# Created: 2019-08-09

# Warning: this script requires sorted dictionary behavior and will likely not
# work with Python < 3.6

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


def profile_reads(fastq_seqs, length):
    position_freq_dict = {}
    for i in range(length):
        position = i + 1
        position_freq_dict.update({position: {'A': 0, 'T': 0, 'C': 0, 'G': 0}})
    for sequence in fastq_seqs:
        position = 0
        if len(sequence) == length:
            for base in sequence:
                position += 1
                position_freq_dict[position][base] += 1
    return position_freq_dict


# Command line parser

def get_args():
    parser = ArgumentParser(
        description='Profile nucleotide content by position from a .fastq file')
    parser.add_argument('fastq',
                        help='Input .fastq, may be gzipped',
                        metavar='FILE.fastq(.gz)')
    parser.add_argument('-l', '--length',
                        help='Length of reads to profile',
                        type=int,
                        metavar='INT')
    return parser.parse_args()


# Main function entry point

def main(args):
    read_profile = profile_reads(fastq_yield_seqs(args.fastq), args.length)
    print('position', 'A', 'T', 'C', 'G', sep=',')
    for position, freq_dict in read_profile.items():
        line = [str(position)]
        for base, count in freq_dict.items():
            line.append(str(count))
        print(','.join(line))


if __name__ == '__main__':
    main(get_args())
