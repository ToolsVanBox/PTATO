process CreateSvPlots {
  tag {"CreateSvPlots ${normal_sample_id} ${tumor_sample_id}"}
  label 'CreateSvPlots'
  shell = ['/bin/bash', '-euo', 'pipefail']

  input:
    tuple( val(donor_id), val(normal_sample_id), val(tumor_sample_id), path(readcounts_100kb_file), path(readcounts_1mb_file), path(readcounts_segments_file), path(baf_binned_100kb_file), path(baf_segments_file), path(cnv_file) )

  output:
    tuple( val(donor_id), val(normal_sample_id), val(tumor_sample_id), path("${tumor_sample_id}.*.p*"), emit: sv_plots)

  script:
    """
    host=\$(hostname)
    echo \${host}

    R --slave --file=${baseDir}/scripts/R/create_sv_plots.R --args \
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
