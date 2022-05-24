process FilterGripss {
  tag {"FilterGripss ${normal_sample_id} ${tumor_sample_id}"}
  label 'FilterGripss'
  shell = ['/bin/bash', '-euo', 'pipefail']

  input:
    tuple( val(donor_id), val(normal_sample_id), val(tumor_sample_id), path(gripss_somatic_filtered_vcf), path(gripss_somatic_filtered_tbi), path(baf_filtered_file), path(cobalt_filtered_readcounts_file) )

  output:
    tuple( val(donor_id), val(normal_sample_id), val(tumor_sample_id), path("${tumor_sample_id}.svs.*"), emit: gripss_filtered_files)

  script:
    """
    host=\$(hostname)
    echo \${host}

    R --slave --file=${baseDir}/scripts/R/filter_gripss.R --args \
    ${gripss_somatic_filtered_vcf} \
    ${baf_filtered_file} \
    ${cobalt_filtered_readcounts_file} \
    ${params.gripss.pon} \
    ${tumor_sample_id} \
    ${params.filtergripss.optional}
    """
}
