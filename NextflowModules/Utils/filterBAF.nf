process FilterBAF {
  tag {"FilterBAF ${normal_sample_id} ${tumor_sample_id}"}
  label 'FilterBAF'
  shell = ['/bin/bash', '-euo', 'pipefail']
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://vanboxtelbioinformatics/ptato_r:1.3.3':
        'europe-west4-docker.pkg.dev/pmc-gcp-box-d-pip-development/pipeline-containers/ptato_r@sha256:c2396d1d9c217123444f4fed846fba38bb5dade14833b244a33bb9806577b059' }"

  input:
    tuple( val(donor_id), val(normal_sample_id), val(tumor_sample_id), path(germline_vcf_file), path(germline_tbi) )
    path( centromeres )
    path( cytoband )

  output:
    tuple( val(donor_id), val(normal_sample_id), val(tumor_sample_id), path("${tumor_sample_id}.baf.*"), emit: baf_filtered_files)

  script:
    """
    host=\$(hostname)
    echo \${host}

    R --slave --file=/scripts/R/filter_baf.R --args \
    ${germline_vcf_file} \
    ${tumor_sample_id} \
    ${normal_sample_id} \
    ${centromeres} \
    ${cytoband} \
    ${params.ref_genome} \
    ${tumor_sample_id} \
    ${params.filterbaf.optional}
    """
}
