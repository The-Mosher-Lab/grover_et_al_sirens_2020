#!/usr/bin/env python3

# Author: Jeffrey Grover
# Purpose: Pull sequences from a fasta file based on coordinates in a bed file
# Created: 2019-06-13

from argparse import ArgumentParser
from itertools import groupby

# Subroutine functions


def fasta_iterate(fasta_file):
    with open(fasta_file, 'r') as input_handle:
        fasta_reader = (
            x[1] for x in groupby(input_handle, lambda line: line.startswith('>'))
        )
        for header in fasta_reader:
            chromosome = next(header).strip('>').strip()
            seq = ''.join(s.strip() for s in next(fasta_reader))
            yield (chromosome, seq)


def parse_bed_to_dict(bed_file):
    bed_dict = {}
    with open(bed_file, 'r') as input_handle:
        for line in input_handle:
            entry = line.strip().split('\t')
            chromosome = entry[0]
            start = int(entry[1])
            stop = int(entry[2]) - 1  # To convert from half-open coordinates
            feature_id = entry[3]
            if chromosome not in bed_dict:
                bed_dict[chromosome] = {}
            bed_dict[chromosome][feature_id] = (start, stop)
    return bed_dict


def wrap_text(text, width=80):
    for s in range(0, len(text), width):
        yield text[s:s+width]


def output_sequences_as_fasta(fasta_iterator, bed_dict):
    for chromosome, sequence in fasta_iterator:
        if chromosome in bed_dict:
            for feature_id in bed_dict[chromosome]:
                start = bed_dict[chromosome][feature_id][0]
                stop = bed_dict[chromosome][feature_id][1]
                header = '>%s:%s-%s_%s' % (chromosome, start + 1, stop + 1, feature_id)
                print(header, sep='\n')
                for line in wrap_text(sequence[start:stop]):
                    print(line)


# CLI argument parser


def get_args():
    parser = ArgumentParser(
        description='Pull sequences from a .fasta file based on coordinates, '
        'using a bed file as input.'
    )
    parser.add_argument('fasta',
                        help='Input .fasta file to use as a reference',
                        metavar='FILE.fasta')
    parser.add_argument('-b', '--bed',
                        help='Input .bed file with at least 4 columns; '
                        'chromosome, start, stop, and ID.',
                        metavar='FILE.bed')
    return parser.parse_args()


# Main function entry point


def main(args):
    bed_dict = parse_bed_to_dict(args.bed)
    output_sequences_as_fasta(fasta_iterate(args.fasta), bed_dict)


if __name__ == '__main__':
    args = get_args()
    main(args)
