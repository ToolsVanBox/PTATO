process createABtable {
  tag {"createABtable ${sample_id} ${chrom}"}
  label 'createABtable'
  shell = ['/bin/bash', '-euo', 'pipefail']
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://vanboxtelbioinformatics/ptato_r:1.3.3':
        'europe-west4-docker.pkg.dev/pmc-gcp-box-d-pip-development/pipeline-containers/ptato_r@sha256:c2396d1d9c217123444f4fed846fba38bb5dade14833b244a33bb9806577b059' }"

  input:
    tuple( val(donor_id), val(sample_id), val(chrom), path(phased_vcf), path(phased_tbi), path(germline_vcf), path(germline_tbi), path(somatic_vcf), path(somatic_tbi) )

  output:
    tuple( val(donor_id), val(sample_id), path("${sample_id}_${chrom}.abtable.txt"), emit: ab_table)

  script:
    """
    host=\$(hostname)
    echo \${host}

    R --slave --file=/scripts/R/ABscript.R --args ${somatic_vcf} ${germline_vcf} ${phased_vcf} ${chrom} ${sample_id}_${chrom}.abtable.txt ${params.ref_genome}
    """
}

process mergeABtable {
  tag {"mergeABtable ${sample_id}"}
  label 'mergeABtable'
  shell = ['/bin/bash', '-euo', 'pipefail']
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://vanboxtelbioinformatics/ptato_r:1.3.3':
        'europe-west4-docker.pkg.dev/pmc-gcp-box-d-pip-development/pipeline-containers/ptato_r@sha256:c2396d1d9c217123444f4fed846fba38bb5dade14833b244a33bb9806577b059' }"

  input:
    tuple( val(donor_id), val(sample_id), path(ab_tables) )

  output:
    tuple( val(donor_id), val(sample_id), path("${sample_id}.abtable.txt"), emit: ab_table )

  script:
    """
    header=true
    for ab_table in ${ab_tables}; do
      if \$header; then
        cat \$ab_table > ${sample_id}.abtable.txt
      else
        cat \$ab_table | tail -n+2 >> ${sample_id}.abtable.txt
      fi
      header=false
    done
    """
}
