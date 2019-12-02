#!/usr/bin/env python3

# Author: Jeffrey Grover
# Purpose: Filter .fastq reads between sizes defined by user input.
# Created: 12/2016

from argparse import ArgumentParser
import gzip


def magic_open(input_file):
    if input_file.endswith('gz'):
        return gzip.open(input_file, 'rt')
    else:
        return open(input_file, 'r')


def filter_by_length(input_path, output_path, min_length, max_length):
    n = 0
    fastq_record = []
    with magic_open(input_path, 'r') as input_file:
        for line in input_file:
            n += 1
            fastq_record.append(line.strip())
            if n == 4:
                if min_length <= len(fastq_record[1]) <= max_length:
                    print('\n'.join(fastq_record))
                n = 0
                fastq_record = []


# Parse command line options

def get_args():
    parser = ArgumentParser(
        description='Filters a given fastq file for reads between a supplied'
        'minimum and maximum length.')
    parser.add_argument('fastq',
                        help='Input .fastq file or fastq.gz',
                        metavar='FILE.fastq(.gz)')
    parser.add_argument('-n', '--min', help='Minimum length for filtering', type=int)
    parser.add_argument('-m', '--max', help='Maximum length for filtering', type=int)
    return parser.parse_args()


# Filter the .fastq

def main(args):
    filter_by_length(args.fastq, args.min, args.max)


if __name__ == '__main__':
    main(get_args())
