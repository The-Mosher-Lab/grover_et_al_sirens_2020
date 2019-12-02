#!/usr/bin/env python3

# Author: Jeffrey Grover
# Purpose: Determine bisulfite conversion rate from MethylDackel bedGraph files
# Created: 2/2019

import gzip
from argparse import ArgumentParser


def magic_open(input_file):
    if input_file.endswith('gz'):
        return gzip.open(input_file, 'rt')
    else:
        return open(input_file, 'r')


def parse_bedgraph(input_context_bedgraph):
    met_count = 0
    unmet_count = 0
    with magic_open(input_context_bedgraph) as input_handle:
        next(input_handle)  # Skip header
        for line in input_handle:
            entry = line.strip().split()
            met_count += int(entry[4])
            unmet_count += int(entry[5])
    return (met_count, unmet_count)


def conversion_calc(cg_counts, chg_counts, chh_counts):
    total_cg = cg_counts[0] + cg_counts[1]
    total_chg = chg_counts[0] + chg_counts[1]
    total_chh = chh_counts[0] + chh_counts[1]
    conversion_rate = (sum(cg_counts[1], chg_counts[1], chh_counts[1]) /
                       sum(total_cg, total_chg, total_chh)) * 100
    return conversion_rate


# Command line parser


def get_args():
    parser = ArgumentParser(
        description=
        'Load CG, CHG, and CHH context bedGraph files containing count of '
        'methylated and unmethylated reads in columns 5 and 6 (ex. from '
        'MethylDackel) and return the conversion rate.')
    parser.add_argument('--CG',
                        help='CG Context bedGraph file.',
                        metavar='FILE.bedGraph(.gz)')
    parser.add_argument('--CHG',
                        help='CHG context bedGraph file.',
                        metavar='FILE.bedGraph(.gz)')
    parser.add_argument('--CHH',
                        help='CHH context bedGraph file.',
                        metavar='FILE.bedGraph(.gz)')
    return parser.parse_args()


# Process the files


def main(args):
    cg_counts = parse_bedgraph(args.CG)
    chg_counts = parse_bedgraph(args.CHG)
    chh_counts = parse_bedgraph(args.CHH)
    conversion_rate = conversion_calc(cg_counts, chg_counts, chh_counts)

    print('CG Methylated/Total:\t', cg_counts[0], '/', sum(cg_counts))
    print('CHG Methylated/Total:\t', chg_counts[0], '/', sum(chg_counts))
    print('CHH Methylated/Total:\t', chh_counts[0], '/', sum(chh_counts))
    print('Conversion Rate:\t', conversion_rate)


if __name__ == '__main__':
    main(get_args())
