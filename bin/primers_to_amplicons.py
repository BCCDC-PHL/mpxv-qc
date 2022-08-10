#!/usr/bin/env python3

import argparse
import csv
import json


def read_primers(primer_bed_path, primer_pairs_path, primer_name_delimiter='_', amplicon_number_offset=2):
    primer_pairs = {}
    with open(primer_pairs_path, 'r') as f:
        for line in f:
            left, right = line.strip().split('\t')
            primer_pairs[left] = right
            primer_pairs[right] = left

    primers_by_name = {}
    with open(primer_bed_path, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            name = fields[3]
            pool = fields[4]
            orientation = fields[5]
            pair_name = primer_pairs[name]
            amplicon_number = name.split(primer_name_delimiter)[amplicon_number_offset]
            primers_by_name[name] = {
                'chrom': chrom,
                'name': name,
                'pair_name': pair_name,
                'start': start,
                'end': end,
                'amplicon_number': amplicon_number,
                'pool': pool,
                'orientation': orientation
            }
            if orientation == '+':
                pass
                
    return primers_by_name


def primers_to_amplicons(primers):
    amplicons = []
    positive_orientation_primers = [primers[p] for p in primers if primers[p]['orientation'] == '+']
    for primer in positive_orientation_primers:
        amplicon_chrom = primer['chrom']
        amplicon_name = '_'.join(primer['pair_name'].split('_')[:-1])
        amplicon_number = int(primer['amplicon_number'])
        amplicon_pool = primer['pool']
        amplicon_start = primer['end']
        amplicon_end = primers[primer['pair_name']]['start']
        amplicon_length = amplicon_end - amplicon_start
        amplicons.append({
            'chrom': amplicon_chrom,
            'name': amplicon_name,
            'number': amplicon_number,
            'start': amplicon_start,
            'end': amplicon_end,
            'length': amplicon_length,
            'pool': amplicon_pool,
        })

    amplicons_sorted = sorted(amplicons, key = lambda x: x['number'])

    return amplicons_sorted

    
def main(args):
    primers = read_primers(args.primer_bed, args.primer_pairs)
    amplicons = primers_to_amplicons(primers)

    for amplicon in amplicons:
        print('\t'.join([
            amplicon['chrom'],
            str(amplicon['start']),
            str(amplicon['end']),
            amplicon['name'],
            amplicon['pool'],
            '+',
        ]))



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--primer-bed')
    parser.add_argument('--primer-pairs')
    args = parser.parse_args()
    main(args)
