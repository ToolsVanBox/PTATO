process IntegrateSvFiles {
  tag {"IntegrateSvFiles ${normal_sample_id} ${tumor_sample_id}"}
  label 'IntegrateSvFiles'
  shell = ['/bin/bash', '-euo', 'pipefail']

  input:
    tuple( val(donor_id), val(normal_sample_id), val(tumor_sample_id), path(baf_filtered_file), path(readcounts_file), path(baf_segments), path(readcounts_segments), path(gripss_filtered_vcf) )

  output:
    tuple( val(donor_id), val(normal_sample_id), val(tumor_sample_id), path("${tumor_sample_id}.integrated.*"), emit: integrated_sv_files)

  script:
    """
    host=\$(hostname)
    echo \${host}

    R --slave --file=${baseDir}/scripts/R/integrate_sv_files.R --args \
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
