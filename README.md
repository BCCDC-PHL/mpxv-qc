# mpxv-qc

A companion pipeline to [BCCDC-PHL/mpxv-artic-nf](https://github.com/BCCDC-PHL/mpxv-artic-nf) for quality-control analysis.
This pipeline is based directly on [jts/ncov-tools](https://github.com/jts/ncov-tools), with some adjustments made for monkeypox virus analysis.

```mermaid
flowchart TD
run_dir[run_dir]
primer_bed[primer.bed]
primer_pairs[primer_pairs.tsv]
run_dir -- consensus --> identify_complete_genomes(identify_complete_genomes)
identify_complete_genomes --> prepare_multi_fasta(prepare_multi_fasta)
nextclade_dataset(nextclade_dataset)
nextclade_dataset -- dataset --> nextclade(nextclade)
prepare_multi_fasta --> augur_align(augur_align)
nextclade_dataset -- ref --> augur_align
augur_align --> augur_tree(augur_tree)
nextclade_dataset -- ref --> augur_tree
augur_align --> make_alleles(make_alleles)
nextclade_dataset -- ref --> make_alleles
augur_tree --> plot_tree_snps(plot_tree_snps)
make_alleles --> plot_tree_snps
nextclade -- qc --> plot_tree_snps
plot_tree_snps --> tree_snps.pdf
run_dir -- variants --> snpeff(snpeff)
nextclade_dataset -- ref --> snpeff
snpeff --> make_aa_table(make_aa_table)
primer_bed --> primer_bed_to_amplicon_bed(primer_bed_to_amplicon_bed)
primer_pairs --> primer_bed_to_amplicon_bed
run_dir -- alignment --> calc_amplicon_depth(calc_amplicon_depth)
primer_bed_to_amplicon_bed --> calc_amplicon_depth
run_dir -- variants --> create_primer_snp_bed(create_primer_snp_bed)
primer_bed --> create_primer_snp_bed
nextclade_dataset -- ref --> make_genome_bed(make_genome_bed)
run_dir -- alignments --> calc_per_base_depth(calc_per_base_depth)
make_genome_bed --> calc_per_base_depth
make_sample_qc_summary(make_sample_qc_summary)
run_dir -- consensus --> make_sample_qc_summary
run_dir -- variants --> make_sample_qc_summary
calc_per_base_depth --> make_sample_qc_summary
make_alleles --> make_sample_qc_summary
make_sample_qc_summary --> summary_qc.csv
```

## Usage
This pipeline makes some assumptions about the directory structure of the input dataset. We assume that there is a directory below the directory passed via the `--run_dir` flag, named like: `mpxv-artic-nf-vX.Y-output`, which contains output from the [BCCDC-PHL/mpxv-artic-nf](https://github.com/BCCDC-PHL/mpxv-artic-nf) pipeline.

```
nextflow run BCCDC-PHL/mpxv-qc \
  -profile conda \
  --cache ~/.conda/envs \
  --run_dir </path/to/run_dir> \
  --bed </path/to/mpxv.scheme.bed> \
  --primer_pairs_tsv </path/to/primer_pairs.tsv> \
  --outdir <output_dir>
```

## Parameters

| Name                          | Default | Description                                                                                      |
|:------------------------------|--------:|:-------------------------------------------------------------------------------------------------|
| `minimum_genome_completeness` | `85.0`  | Genome completeness threshold below which samples will be excluded from tree/SNPs plot.          |
| `partial_genome_threshold`    | `85.0`  | Genome completeness threshold below which samples will be tagged as `PARITAL_GENOME`             |
| `incomplete_genome_threshold` | `50.0`  | Genome completeness threshold below which samples will be tagged as `INCOMPLETE_GENOME`          |
| `excess_ambiguity_threshold`  | `5`     | Number of ambiguous bases in consensus above which samples will be tagged as `EXCESS_AMBIGUITY`  |

