process FilterCobalt {
  tag {"FilterCobalt ${normal_sample_id} ${tumor_sample_id}"}
  label 'FilterCobalt'
  shell = ['/bin/bash', '-euo', 'pipefail']

  input:
    tuple( val(donor_id), val(normal_sample_id), val(tumor_sample_id), path(cobalt_tsv_file) )

  output:
    tuple( val(donor_id), val(normal_sample_id), val(tumor_sample_id), path("${tumor_sample_id}.readcounts.*"), emit: cobalt_filtered_files)

  script:
    """
    host=\$(hostname)
    echo \${host}

    R --slave --file=${baseDir}/scripts/R/filter_cobalt.R --args \
    ${cobalt_tsv_file} \
    ${params.svs.centromeres} \
    ${params.svs.cytoband} \
    ${params.cobalt.pon} \
    ${params.ref_genome} \
    ${tumor_sample_id} \
    ${params.filtercobalt.optional}
    """
}
