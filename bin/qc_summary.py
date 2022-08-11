#!/usr/bin/env python3

import argparse
import csv
import statistics
import sys

import pysam


def count_total_bases_in_fasta(fasta_path):
    fasta = pysam.FastxFile(fasta_path)
    sequence = ""

    for record in fasta:
        sequence = record.sequence.upper()

    total_bases = 0

    for base in sequence:
        total_bases += 1

    return total_bases


def count_bases(sequence, bases):
    num_bases = 0

    for base in sequence:
        if base in bases:
            num_bases += 1

    return num_bases


def get_variant_counts(vcf_path):
    vcf = pysam.VariantFile(vcf_path, 'r')

    num_snvs = 0
    num_indel = 0
    num_indel_triplet = 0

    for rec in vcf.fetch():
        record_type = None
        record_alt = rec.alts[0]
        record_ref_len = len(rec.ref)
        record_alt_len = len(record_alt)
        indel_len = abs(record_ref_len - record_alt_len)

        for k, v in rec.info.items():
            if k == 'TYPE':
                if v[0] == 'snp':
                    record_type = 'snp'
                elif v[0] == 'ins' or v[0] == 'del':
                    record_type = 'indel'

        if record_type == 'snp' and indel_len == 0:
            num_snvs += 1

        if record_type == 'indel':
            num_indel += 1

        if record_type == 'indel' and indel_len % 3 == 0:
            num_indel_triplet += 1

    variant_counts = {
        'num_snvs': num_snvs,
        'num_indel': num_indel,
        'num_indel_triplet': num_indel_triplet,
    }

    return variant_counts



def get_num_consensus_snvs(alleles_tsv_path, sample_id, excluded_alleles):

    num_consensus_snvs = 0
    
    with open(alleles_tsv_path, 'r') as f:
        reader = csv.DictReader(f, dialect='excel-tab')
        for row in reader:
            if row['name'] == sample_id:
                if row['alt_allele'] not in excluded_alleles:
                    num_consensus_snvs += 1

    return num_consensus_snvs


def get_coverage_stats(per_base_coverage_path, delimiter='\t'):
    coverage_stats = {}
    depth = []
    with open(per_base_coverage_path) as f:
        reader = csv.DictReader(f, delimiter=delimiter)
        for row in reader:
            depth.append(int(row['depth']))

    mean_depth = round(statistics.mean(depth), 1)
    median_depth = round(statistics.median(depth), 1)

    coverage_stats['mean_sequencing_depth'] = mean_depth
    coverage_stats['median_sequencing_depth'] = median_depth
    
    return coverage_stats


def main(args):

    iupac_codes = ['R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V']

    qc_line = {}
    qc_line['sample'] = args.sample
    qc_line['run_name'] = args.run_name

    consensus_fasta = pysam.FastxFile(args.consensus)
    
    consensus_seq = ""
    # There should be only one record in the consensus fasta
    for record in consensus_fasta:
        consensus_seq = record.sequence.upper()

    total_consensus_bases = len(consensus_seq)

    qc_line['num_consensus_snvs'] = get_num_consensus_snvs(args.alleles, args.sample, iupac_codes + ['N'])

    num_consensus_n = count_bases(consensus_seq, ['N'])
    
    qc_line['num_consensus_n'] = num_consensus_n

    if total_consensus_bases > 0:
        qc_line['genome_completeness'] = 1 - round(num_consensus_n / float(total_consensus_bases), 4)

    qc_line['num_consensus_iupac'] = count_bases(consensus_seq, iupac_codes)

    variant_counts = get_variant_counts(args.variants)

    qc_line['num_variants_snvs'] = variant_counts['num_snvs']

    qc_line['num_variants_indel'] = variant_counts['num_indel']

    qc_line['num_variants_indel_triplet'] = variant_counts['num_indel_triplet']
    
    coverage_stats = get_coverage_stats(args.coverage)

    qc_line['mean_sequencing_depth'] = coverage_stats['mean_sequencing_depth']

    qc_line['median_sequencing_depth'] = coverage_stats['median_sequencing_depth']
    
    output_fieldnames = [
        'sample',
        'run_name',
        'num_consensus_snvs',
        'num_consensus_n',
        'num_consensus_iupac',
        'num_variants_snvs',
        'num_variants_indel',
        'num_variants_indel_triplet',
        'mean_sequencing_depth',
        'median_sequencing_depth',
        'genome_completeness',
        'qc_pass',
    ]
    
    writer = csv.DictWriter(sys.stdout, fieldnames=output_fieldnames, dialect='excel-tab')
    writer.writeheader()
    writer.writerow(qc_line)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Tool for summarizing QC data")
    parser.add_argument('-c', '--consensus', help='<sample>.consensus.fasta file to process')
    parser.add_argument('-v', '--variants',
                        help='<sample>.vcf file to process')
    parser.add_argument('-e', '--coverage',
                        help='<sample>.per_base_coverage.bed file to process')
    parser.add_argument('-i', '--indel', action='store_true',
                        help='flag to determine whether to count indels')
    parser.add_argument('-m', '--meta', default=None,
                        help='full path to the metadata file')
    parser.add_argument('-a', '--alleles',
                        help='full path to the alleles.tsv file')
    parser.add_argument('-s', '--sample',
                        help='name of sample being processed')
    parser.add_argument('-x', '--mixture', default=None,
                        help='full path to the mixture report')
    parser.add_argument('-r', '--run_name',
                        help='run name for sample')
    parser.add_argument('-t', '--aa_table',
                        help='full path to the <sample>_aa_table.tsv file')

    args = parser.parse_args()

    main(args)
