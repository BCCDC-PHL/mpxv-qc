#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { nextclade_dataset } from './modules/mpxv-qc.nf'
include { identify_complete_genomes } from './modules/mpxv-qc.nf'
include { prepare_multi_fasta } from './modules/mpxv-qc.nf'
include { nextclade } from './modules/mpxv-qc.nf'
include { augur_align } from './modules/mpxv-qc.nf'
include { augur_tree } from './modules/mpxv-qc.nf'
include { make_alleles } from './modules/mpxv-qc.nf'
include { plot_tree_snps } from './modules/mpxv-qc.nf'

workflow {

  ch_run_dir            = Channel.fromPath(params.run_dir)

  ch_run_id             = Channel.of(file(params.run_dir).getName())

  ch_artic_analysis_dir = Channel.fromPath(params.run_dir + "/" + params.artic_analysis_subdir)

  main:

    identify_complete_genomes(ch_run_id.combine(ch_artic_analysis_dir))

    prepare_multi_fasta(identify_complete_genomes.out.combine(ch_artic_analysis_dir))

    nextclade_dataset(ch_run_id)

    nextclade(prepare_multi_fasta.out.join(nextclade_dataset.out.dataset))

    augur_align(prepare_multi_fasta.out.join(nextclade_dataset.out.ref))

    augur_tree(augur_align.out.join(nextclade_dataset.out.ref))

    make_alleles(augur_align.out.join(nextclade_dataset.out.ref))

    plot_tree_snps(augur_tree.out.join(make_alleles.out).join(nextclade.out.qc))
     
}
