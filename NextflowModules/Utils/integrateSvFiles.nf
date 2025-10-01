process IntegrateSvFiles {
  tag {"IntegrateSvFiles ${normal_sample_id} ${tumor_sample_id}"}
  label 'IntegrateSvFiles'
  shell = ['/bin/bash', '-euo', 'pipefail']
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?  
        'docker://vanboxtelbioinformatics/ptato_r:1.3.3':
        'europe-west4-docker.pkg.dev/pmc-gcp-box-d-pip-development/pipeline-containers/ptato_r@sha256:c2396d1d9c217123444f4fed846fba38bb5dade14833b244a33bb9806577b059' }"

  input:
    tuple( val(donor_id), val(normal_sample_id), val(tumor_sample_id), path(baf_filtered_file), path(readcounts_file), path(baf_segments), path(readcounts_segments), path(gripss_filtered_vcf) )

  output:
    tuple( val(donor_id), val(normal_sample_id), val(tumor_sample_id), path("${tumor_sample_id}.integrated.*"), emit: integrated_sv_files)

  script:
    """
    host=\$(hostname)
    echo \${host}

    R --slave --file=/scripts/R/integrate_sv_files.R --args \
    ${baf_filtered_file} \
    ${readcounts_file} \
    ${baf_segments} \
    ${readcounts_segments} \
    ${gripss_filtered_vcf} \
    ${tumor_sample_id} \
    ${tumor_sample_id} \
    ${params.integratesvfiles.optional}
    """
}
