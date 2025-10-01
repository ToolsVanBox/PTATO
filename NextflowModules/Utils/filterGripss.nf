process FilterGripss {
  tag {"FilterGripss ${normal_sample_id} ${tumor_sample_id}"}
  label 'FilterGripss'
  shell = ['/bin/bash', '-euo', 'pipefail']
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?  
        'docker://vanboxtelbioinformatics/ptato_r:1.3.3':
        'europe-west4-docker.pkg.dev/pmc-gcp-box-d-pip-development/pipeline-containers/ptato_r@sha256:c2396d1d9c217123444f4fed846fba38bb5dade14833b244a33bb9806577b059' }"

  input:
    tuple( val(donor_id), val(normal_sample_id), val(tumor_sample_id), path(gripss_somatic_filtered_vcf), path(gripss_somatic_filtered_tbi), path(baf_filtered_file), path(cobalt_filtered_readcounts_file) )
    path( pon )
  output:
    tuple( val(donor_id), val(normal_sample_id), val(tumor_sample_id), path("${tumor_sample_id}.svs.*"), emit: gripss_filtered_files)

  script:
    """
    host=\$(hostname)
    echo \${host}

    R --slave --file=/scripts/R/filter_gripss.R --args \
    ${gripss_somatic_filtered_vcf} \
    ${baf_filtered_file} \
    ${cobalt_filtered_readcounts_file} \
    ${pon} \
    ${tumor_sample_id} \
    ${params.filtergripss.optional}
    """
}
