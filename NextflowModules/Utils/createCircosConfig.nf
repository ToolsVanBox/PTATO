process CreateCircosConfig {
  tag {"CreateCircosConfig ${normal_sample_id} ${tumor_sample_id}"}
  label 'CreateCircosConfig'
  shell = ['/bin/bash', '-euo', 'pipefail']
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://vanboxtelbioinformatics/ptato_r:1.3.3':
        'europe-west4-docker.pkg.dev/pmc-gcp-box-d-pip-development/pipeline-containers/ptato_rsha256:01608f437f46324b7ef4c6ae3ce9fe06da8d160ab6d46a44b04342e979f13be6' }"

  input:
    tuple( val(donor_id), val(normal_sample_id), val(tumor_sample_id), path(cnv_file), path(readcounts_1mb_file), path(baf_1mb_file), path(readcounts_segments_txt_file), path(baf_segments_txt_file), path(sv_vcf) )

  output:
    tuple( val(donor_id), val(normal_sample_id), val(tumor_sample_id), path("${tumor_sample_id}.circos*"), emit: circos_config_files)

  script:
    """
    host=\$(hostname)
    echo \${host}

    R --slave --file=/scripts/R/create_circos_config.R --args \
    ${cnv_file} \
    ${readcounts_1mb_file} \
    ${baf_1mb_file} \
    ${readcounts_segments_txt_file} \
    ${baf_segments_txt_file} \
    ${sv_vcf} \
    ${tumor_sample_id} \
    ${params.createcircosconfig.optional}
    """
}
