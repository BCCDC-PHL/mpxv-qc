#!/usr/bin/env python3

import argparse
import csv
import json
import re
import sys

import pysam

relevant_annotations = [
    'missense_variant',
    'synonymous_variant',
    'frameshift_variant',
    'start_lost',
    'stop_gained',
    'conservative_inframe_deletion',
    'conservative_inframe_insertion',
    'disruptive_inframe_deletion',
    'disruptive_inframe_insertion',
]


def process_variant_record(var, sample, relevant_annotations):
    var_meta = {}
    
    for k, v in var.info.items():
        if k == 'ANN':
            for annotations in v:
                annotation = annotations.split('|')
                if annotation[1] in relevant_annotations:
                    var_meta['sample'] = sample
                    var_meta['chr'] = var.chrom
                    var_meta['pos'] = var.pos
                    var_meta['ref'] = var.ref
                    var_meta['alt'] = ','.join(var.alts)
                    var_meta['Consequence'] = annotation[1]
                    var_meta['gene'] = annotation[3]
                    var_meta['protein'] = re.sub('^p.', '', annotation[10])

                    if var_meta['gene'] != '' and var_meta['protein'] != '':
                        var_meta['aa'] = '-'.join([var_meta['gene'], var_meta['protein']])
                    else:
                        var_meta['aa'] = 'NA'

    return var_meta


def main(args):
    vcf = pysam.VariantFile(args.vcf, 'r')

    vars = []
    for variant_record in vcf.fetch():
        var = process_variant_record(variant_record, args.sample, relevant_annotations)
        vars.append(var)
        
    output_fieldnames = [
        'sample',
        'chr',
        'pos',
        'ref',
        'alt',
        'Consequence',
        'gene',
        'protein',
        'aa'
    ]

    writer = csv.DictWriter(sys.stdout, fieldnames=output_fieldnames, dialect='excel-tab')
    writer.writeheader()
    for var in vars:
        if var:
            writer.writerow(var)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--vcf',
                        help='path to the annotated VCF file from SNPEff')
    parser.add_argument('-s', '--sample',
                        help='name of the sample to be processed')
    args = parser.parse_args()
    main(args)

