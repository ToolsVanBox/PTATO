process CreateSvPlots {
  tag {"CreateSvPlots ${normal_sample_id} ${tumor_sample_id}"}
  label 'CreateSvPlots'
  shell = ['/bin/bash', '-euo', 'pipefail']
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://vanboxtelbioinformatics/ptato_r:1.3.3':
        params.artifact_registry_path + '/ptato_r:1.3.3' }"
//        'europe-west4-docker.pkg.dev/pmc-gcp-box-d-pip-development/pipeline-containers/ptato_r@sha256:c2396d1d9c217123444f4fed846fba38bb5dade14833b244a33bb9806577b059' }"

  input:
    tuple( val(donor_id), val(normal_sample_id), val(tumor_sample_id), path(readcounts_100kb_file), path(readcounts_1mb_file), path(readcounts_segments_file), path(baf_binned_100kb_file), path(baf_segments_file), path(cnv_file) )

  output:
    tuple( val(donor_id), val(normal_sample_id), val(tumor_sample_id), path("${tumor_sample_id}.*.p*"), emit: sv_plots)

  script:
    """
    host=\$(hostname)
    echo \${host}

    R --slave --file=/scripts/R/create_sv_plots.R --args \
    ${readcounts_100kb_file} \
    ${readcounts_1mb_file} \
    ${readcounts_segments_file} \
    ${baf_binned_100kb_file} \
    ${baf_segments_file} \
    ${cnv_file} \
    ${tumor_sample_id} \
    ${tumor_sample_id} \
    ${params.createsvplots.optional}
    """
}
