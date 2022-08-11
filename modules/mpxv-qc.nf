process identify_complete_genomes {

  tag { run_id }

  executor 'local'

  input:
    tuple val(run_id), path(artic_analysis_dir)

  output:
    tuple val(run_id), path("${run_id}_qc_pass_sample_ids.csv")

  script:
  """
  tail -qn+2 ${artic_analysis_dir}/${run_id}.qc.csv | awk -F ',' '\$3 > (${params.minimum_genome_completeness} * 100)' | cut -d ',' -f 1 | grep -v '^POS' | grep -v '^NEG' > ${run_id}_qc_pass_sample_ids.csv
  """
}


process prepare_multi_fasta {

  tag { run_id }

  executor 'local'

  input:
    tuple val(run_id), path(qc_pass_sample_ids), path(artic_analysis_dir)

  output:
    tuple val(run_id), path("${run_id}_qc_pass_seqs.fa")

  script:
  """
  while read -r sample_id; do
    cat ${artic_analysis_dir}/${params.artic_consensus_subdir}/\${sample_id}.consensus.fa >> ${run_id}_qc_pass_seqs.fa
  done < ${qc_pass_sample_ids}
  """
}


process nextclade_dataset {

  tag { run_id }

  executor 'local'

  input:
    val(run_id)

  output:
    tuple val(run_id), path("nextclade_${params.nextclade_dataset}"), emit: dataset
    tuple val(run_id), path("nextclade_${params.nextclade_dataset}_ref.fasta"), emit: ref
    tuple val(run_id), path("nextclade_${params.nextclade_dataset}_genemap.gff"), emit: genemap

  script:
  """
  nextclade dataset get --name ${params.nextclade_dataset} --output-dir nextclade_${params.nextclade_dataset}
  cp nextclade_${params.nextclade_dataset}/reference.fasta nextclade_${params.nextclade_dataset}_ref.fasta
  cp nextclade_${params.nextclade_dataset}/genemap.gff nextclade_${params.nextclade_dataset}_genemap.gff
  """
}


process nextclade {

  tag { run_id }

  publishDir "${params.outdir}", pattern: "${run_id}_nextclade_qc.tsv", mode: 'copy'

  input:
    tuple val(run_id), path(sequences), path(dataset)

  output:
    tuple val(run_id), path("${run_id}.aln.fa"), emit: alignment
    tuple val(run_id), path("${run_id}_nextclade_qc.tsv"), emit: qc

  script:
  """
  nextclade run \
    --jobs ${task.cpus} \
    --input-dataset ${dataset} \
    --include-reference \
    --output-all ${run_id}_nextclade \
    --output-tsv ${run_id}_nextclade_qc.tsv \
    --output-fasta ${run_id}.aln.fa \
    --output-basename ${run_id} \
    ${sequences} \
    > nextclade.log 2>&1
  """
}

process augur_align {

  tag { run_id }

  input:
    tuple val(run_id), path(sequences), path(ref)

  output:
    tuple val(run_id), path("${run_id}_aligned.fasta")

  script:
  """
  augur align \
    --fill-gaps \
    --sequences ${sequences} \
    --reference-sequence ${ref} \
    --output ${run_id}_aligned.fasta
  """
}

process augur_tree {

  tag { run_id }

  publishDir "${params.outdir}/qc_analysis", pattern: "${run_id}_tree.nwk", mode: 'copy'

  input:
    tuple val(run_id), path(alignment), path(ref)

  output:
    tuple val(run_id), path("${run_id}_tree.nwk")

  script:
  """
  augur tree \
    --alignment ${alignment} \
    --output ${run_id}_tree_raw.nwk
  nw_reroot ${run_id}_tree_raw.nwk `head -1 ${ref} | tr -d \">\" | cut -f 1` > ${run_id}_tree.nwk
  """
}

process make_alleles {

  tag { run_id }

  publishDir "${params.outdir}/qc_analysis", pattern: "${run_id}_alleles.tsv", mode: 'copy'

  input:
    tuple val(run_id), path(alignment), path(ref)

  output:
    tuple val(run_id), path("${run_id}_alleles.tsv")


  script:
  """
  align2alleles.py --reference-name `head -1 ${ref} | tr -d \">\" | cut -f 1 -d \" \"` ${alignment} > ${run_id}_alleles.tsv
  """
}

process plot_tree_snps {

  tag { run_id }

  publishDir "${params.outdir}/plots", pattern: "${run_id}_tree_snps.pdf", mode: 'copy'

  input:
    tuple val(run_id), path(tree), path(alleles), path(nextclade_qc)

  output:
    tuple val(run_id), path("${run_id}_tree_snps.pdf")

  script:
  """
  plot_tree_snps.R \
    ${run_id}_tree_snps.pdf \
    ${tree} \
    ${alleles} \
    ${nextclade_qc}
  """
}

process build_snpeff_db {

  input:
    path(ref)

  output:

  script:
  """
  build_snpeff_db.py --accession `head -1 ${ref} | tr -d \">\" | cut -f 1 -d \" \"`
  """
}

process snpeff {

  tag { sample_id }

  publishDir "${params.outdir}/qc_annotation", pattern: "${sample_id}.ann.vcf", mode: 'copy'

  input:
    tuple val(sample_id), path(variants), path(ref)

  output:
    tuple val(sample_id), path("${sample_id}.ann.vcf")

  script:
  """
  snpEff -noLog -hgvs1LetterAa `head -1 ${ref} | tr -d \">\" | cut -f 1 -d \" \"` ${variants} > ${sample_id}.ann.vcf
  """
}

process make_aa_table {

  tag { sample_id }

  publishDir "${params.outdir}/qc_annotation", pattern: "${sample_id}_aa_table.tsv", mode: 'copy'

  input:
    tuple val(sample_id), path(annotated_variants)

  output:
    tuple val(sample_id), path("${sample_id}_aa_table.tsv")

  script:
  """
  snpeff_vcf_to_aa_table.py -s ${sample_id} --vcf ${annotated_variants} > ${sample_id}_aa_table.tsv
  """
}

process primer_bed_to_amplicon_bed {

  tag { primer_bed.baseName }

  publishDir "${params.outdir}/bed", pattern: "amplicon.bed", mode: 'copy'

  input:
    tuple path(primer_bed), path(primer_pairs)

  output:
    path("amplicon.bed")

  script:
  """
  primers_to_amplicons.py --primer-bed ${primer_bed} --primer-pairs ${primer_pairs} > amplicon.bed
  """
}

process make_genome_bed {

  tag { ref.baseName }

  publishDir "${params.outdir}/bed", pattern: "genome.bed", mode: 'copy'

  input:
    tuple val(run_id), path(ref)

  output:
    path("genome.bed")

  script:
  """
  samtools faidx ${ref}
  cat ${ref}.fai | awk '{{ print \$1 \"\\t0\\t\" \$2 }}' > genome.bed
  """
}


process calc_amplicon_depth {

  tag { sample_id }

  publishDir "${params.outdir}/qc_sequencing", pattern: "${sample_id}.amplicon_depth.bed", mode: 'copy'

  input:
    tuple val(sample_id), path(alignment), path(alignment_index), path(amplicon_bed)

  output:
    tuple val(sample_id), path("${sample_id}.amplicon_depth.bed")

  script:
  """
  echo -e \"reference_name\tstart\tend\tamplicon_id\tpool\tstrand\tmean_depth\" > ${sample_id}.amplicon_depth.bed
  bedtools coverage -mean -a ${amplicon_bed} -b ${alignment} >> ${sample_id}.amplicon_depth.bed
  """
}

process calc_per_base_depth {

  tag { sample_id }

  publishDir "${params.outdir}/qc_sequencing", pattern: "${sample_id}.per_base_coverage.bed", mode: 'copy'

  input:
    tuple val(sample_id), path(alignment), path(alignment_index), path(genome_bed)

  output:
    tuple val(sample_id), path("${sample_id}.per_base_coverage.bed")

  script:
  """
  echo -e \"reference_name\tstart\tend\tposition\tdepth" > ${sample_id}.per_base_coverage.bed
  bedtools coverage -d -a ${genome_bed} -b ${alignment} >> ${sample_id}.per_base_coverage.bed
  """
}

process create_primer_snp_bed {

  tag { sample_id }

  publishDir "${params.outdir}/qc_annotation", pattern: "${sample_id}.primer_snps.bed", mode: 'copy'

  input:
    tuple val(sample_id), path(variants), path(primer_bed)

  output:
    tuple val(sample_id), path("${sample_id}.primer_snps.bed")

  script:
  """
  bedtools intersect -a ${variants} -b ${primer_bed} -wb > ${sample_id}.primer_snps.bed
  """
}


process make_sample_qc_summary {

  tag { sample_id }

  input:
    tuple val(sample_id), path(consensus), path(variants), path(per_base_coverage_bed), path(aa_table), val(run_id), path(alleles_tsv)

  output:
    tuple val(sample_id), path("${sample_id}_summary_qc.tsv")

  script:
  """
  qc_summary.py \
    --sample ${sample_id} \
    --run_name ${run_id} \
    --consensus ${consensus} \
    --alleles ${alleles_tsv} \
    --variants ${variants} \
    --coverage ${per_base_coverage_bed} \
    --aa_table ${aa_table} \
    --incomplete_genome_threshold ${params.incomplete_genome_threshold} \
    --partial_genome_threshold ${params.partial_genome_threshold} \
    --excess_ambiguity_threshold ${params.excess_ambiguity_threshold} \
    > ${sample_id}_summary_qc.tsv
  """
}

process write_qc_summary {

  tag { run_id }

  executor 'local'

  input:
    tuple val(run_id), val(lines)

  exec:
  file("${params.outdir}/${run_id}_summary_qc.tsv").withWriter { writer ->
      for ( line in lines ) {
          writer.writeLine(line.join('\t'))
         }   
  }
}