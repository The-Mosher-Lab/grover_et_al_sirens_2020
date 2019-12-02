#!/usr/bin/env python3

# Author: Jeffrey Grover
# Purpose: Calculate percent methylation per feature in a .bed file
# Created: 2019-07-12
# Depends: GNU Sort, bedtools

import subprocess
from argparse import ArgumentParser
from sys import stderr
from os import getpid, remove


# Sorting and calculation functions

def sort_bedgraph(input_file, output_file):
    subprocess.run(
        "sed '1d' {} | sort -k1,1 -k2,2n > {}".format(input_file, output_file),
        shell=True
    )


def sort_features_bed(input_file, output_file):
    subprocess.run(
        'sort -k1,1 -k2,2n {} > {}'.format(input_file, output_file),
        shell=True
    )


def bedtools_map_sum(input_bed, input_bedgraph, output_file):
    subprocess.run(
        'bedtools map -a {} -b {} -c 5,6 -o sum,sum -null 0 > {}'.format(
            input_bed, input_bedgraph, output_file
        ), shell=True
    )


def calc_methylation(input_file, mincov):
    with open(input_file, 'r') as input_handle:
        for line in input_handle:
            entry = line.strip().split()
            nC = int(entry[4])
            nT = int(entry[5])
            if (nC + nT) >= mincov:
                try:
                    perc_met = nC / (nC + nT) * 100
                except ZeroDivisionError:
                    perc_met = 'NA'  # Indicate missing data for zero depth
                print(*entry[0:4], perc_met, sep='\t')


# Get command line options

def get_args():
    parser = ArgumentParser(
        description='Calculate percent methylation over features of interest '
        'from an input .bed file of regions and the per-base or region '
        'methylation calls from a .bedGraph file (from MethylDackel for '
        'example).')
    parser.add_argument('-b', '--bed',
                        help='bed file to process, three columns only',
                        metavar='FILE.bed')
    parser.add_argument('-g', '--bedGraph',
                        help='bedGraph file of methylation calls',
                        metavar='FILE.bedGraph')
    parser.add_argument('-m', '--mincov',
                        help='Minimum coverage value to report methylation',
                        metavar='INT',
                        type=int,
                        default=0)
    parser.add_argument('-s', '--sorted',
                        help='Input files have been pre-sorted',
                        action='store_true')
    parser.add_argument('-k', '--keep_sorted',
                        help='Keep the sorted intermediate files',
                        action='store_true')
    return parser.parse_args()


# Run

def main(args):
    if not args.sorted:
        print('Sorting MethylDackel .bedGraph file: %s' % args.bedGraph, file=stderr)
        sorted_bedGraph = args.bedGraph.replace('.bedGraph', '.sorted.bedGraph')
        sort_bedgraph(args.bedGraph, sorted_bedGraph)

        print('Sorting features .bed file: %s' % args.bed, file=stderr)
        sorted_features_bed = args.bed.replace('.bed', '.sorted.bed')
        sort_features_bed(args.bed, sorted_features_bed)
    else:
        sorted_bedGraph = args.bedGraph
        sorted_features_bed = args.bed

    print('Summing methylation over features...', file=stderr)
    sum_filename = 'bedtools_map_sum.%s.tmp' % getpid()
    bedtools_map_sum(sorted_features_bed, sorted_bedGraph, sum_filename)

    print('Calculating percent methylation per feature...', file=stderr)
    calc_methylation(sum_filename, args.mincov)

    print('Cleaning up temporary files...', file=stderr)
    if not args.keep_sorted and not args.sorted:
        remove(sorted_features_bed)
        remove(sorted_bedGraph)
    remove(sum_filename)


if __name__ == '__main__':
    main(get_args())
