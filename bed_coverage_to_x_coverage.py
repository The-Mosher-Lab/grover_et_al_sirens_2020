#!/usr/bin/env python3

# Author: Jeffrey Grover
# Purpose: Calculate the average coverage over an entire genome from mosdepth
# per-base coverage output
# Created: 2019-07-10

import gzip
from argparse import ArgumentParser


def magic_open(input_file):
    if input_file.endswith('gz'):
        return gzip.open(input_file, 'rt')
    else:
        return open(input_file, 'r')


def get_genome_size(fasta_file):
    genome_size = 0
    with magic_open(fasta_file) as fasta_reader:
        for line in fasta_reader:
            if not line.startswith('>'):
                genome_size += len(line.strip())
    return genome_size


def get_depth(coverage_bedgz):
    bases_sequenced = 0
    with magic_open(coverage_bedgz) as bed_reader:
        for line in bed_reader:
            length = int(line.split('\t')[2]) - int(line.split('\t')[1])
            bases_sequenced += int(line.split('\t')[3]) * length
    return bases_sequenced


def get_x_coverage(bases_sequenced, genome_size):
    return bases_sequenced / genome_size


# Command line Parser

def get_args():
    parser = ArgumentParser(
        description='Calculate X coverage from a bed file with per-base or '
        'window coverage. As output by mosdepth, for example.')
    parser.add_argument(
        '-f', '--fasta',
        help='.fasta file for genome',
        metavar='FILE.fasta')
    parser.add_argument(
        '-m', '--mosdepth',
        help='.bed.gz mosdepth output',
        metavar='FILE.bed.gz')
    return parser.parse_args()


# Process the files

def main(args):
    genome_size = get_genome_size(args.fasta)
    depth = get_depth(args.mosdepth)
    x_coverage = get_x_coverage(depth, genome_size)

    # Output to stdout

    print('Genome Size:', genome_size)
    print('Total Depth:', depth)
    print('X Coverage:', x_coverage)


if __name__ == '__main__':
    main(get_args())
