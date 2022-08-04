process identify_complete_genomes {

  tag { run_id }

  executor 'local'

  input:
    tuple val(run_id), path(artic_analysis_dir)

  output:
    tuple val(run_id), path("${run_id}_qc_pass_sample_ids.csv")

  script:
  """
  tail -qn+2 ${artic_analysis_dir}/${run_id}.qc.csv | awk -F ',' '\$3 > ${params.minimum_genome_completeness}' | cut -d ',' -f 1 | grep -v '^POS' | grep -v '^NEG' > ${run_id}_qc_pass_sample_ids.csv
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

  script:
  """
  nextclade dataset get --name ${params.nextclade_dataset} --output-dir nextclade_${params.nextclade_dataset}
  cp nextclade_${params.nextclade_dataset}/reference.fasta nextclade_${params.nextclade_dataset}_ref.fasta
  """
}


process nextclade {

  tag { run_id }

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

