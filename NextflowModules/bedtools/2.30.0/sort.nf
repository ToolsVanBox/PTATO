process sort {
  tag {"bedtoolsSort ${sample_id}"}
  label 'bedtoolsSort'
  label 'bedtools_2_30_0_sort'
  shell = ['/bin/bash', '-euo', 'pipefail']
  container = 'quay.io/biocontainers/bedtools:2.30.0--h468198e_3'

  input:
    tuple( val(donor_id), val(sample_id), path(bed))

  output:
    tuple( val(donor_id), val(sample_id), path("${sample_id}.sorted.bed"), emit: sorted_bed)

  script:
    """
    host=\$(hostname)
    echo \${host}

    bedtools \
    sort \
    -i ${bed} \
    ${params.bedtoolssort.optional} \
    > ${sample_id}.sorted.bed
    """
}
