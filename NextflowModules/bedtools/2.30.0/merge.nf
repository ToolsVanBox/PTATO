process merge {
  tag {"bedtoolsMerge ${sample_id} ${feature_id}"}
  label 'bedtoolsMerge'
  label 'bedtools_2_30_0_merge'
  shell = ['/bin/bash', '-euo', 'pipefail']
  container = 'quay.io/biocontainers/bedtools:2.30.0--h468198e_3'

  input:
    tuple( val(donor_id), val(sample_id), val(feature_id), path(bed), val(merge_params) )

  output:
    tuple( val(donor_id), val(sample_id), path("${sample_id}.${feature_id}.merged.bed"), emit: feature_merged_bed)

  script:
    """
    host=\$(hostname)
    echo \${host}

    bedtools \
    merge \
    -i ${bed} \
    ${merge_params} \
    ${params.bedtoolsmerge.optional} \
    > ${sample_id}.${feature_id}.merged.bed
    """
}

process mergeAll {
  tag {"bedtoolsMergeAll ${sample_id}"}
  label 'bedtoolsMergeAll'
  label 'bedtools_2_30_0_mergeAll'
  shell = ['/bin/bash', '-euo', 'pipefail']
  container = 'quay.io/biocontainers/bedtools:2.30.0--h468198e_3'

  input:
    tuple( val(donor_id), val(sample_id), path(bed))

  output:
    tuple( val(donor_id), val(sample_id), path("${sample_id}.features.merged.bed"), emit: feature_merged_bed)

  script:
    """
    host=\$(hostname)
    echo \${host}

    bedtools \
    merge \
    -i ${bed} \
    ${params.bedtoolsmergeall.optional} \
    > ${sample_id}.features.merged.bed
    """
}
