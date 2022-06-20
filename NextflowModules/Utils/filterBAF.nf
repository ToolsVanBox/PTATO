process FilterBAF {
  tag {"FilterBAF ${normal_sample_id} ${tumor_sample_id}"}
  label 'FilterBAF'
  shell = ['/bin/bash', '-euo', 'pipefail']

  input:
    tuple( val(donor_id), val(normal_sample_id), val(tumor_sample_id), path(germline_vcf_file), path(germline_tbi) )

  output:
    tuple( val(donor_id), val(normal_sample_id), val(tumor_sample_id), path("${tumor_sample_id}.baf.*"), emit: baf_filtered_files)

  script:
    """
    host=\$(hostname)
    echo \${host}

    R --slave --file=${baseDir}/scripts/R/filter_baf.R --args \
    ${germline_vcf_file} \
    ${tumor_sample_id} \
    ${normal_sample_id} \
    ${params.svs.centromeres} \
    ${params.svs.cytoband} \
    ${params.ref_genome} \
    ${tumor_sample_id} \
    ${params.filterbaf.optional}
    """
}
