process groupby {
  tag {"bedtoolsGroupby ${sample_id} ${feature_id}"}
  label 'bedtoolsGroupby'
  label 'bedtools_2_30_0_groupby'
  shell = ['/bin/bash', '-euo', 'pipefail']
  container = 'quay.io/biocontainers/bedtools:2.30.0--h468198e_3'

  input:
    tuple( val(donor_id), val(sample_id), val(feature_id), path(bed), val(groupby_params) )

  output:
    tuple( val(donor_id), val(sample_id), path("${sample_id}.${feature_id}.groupby.bed"), emit: feature_groupby_bed)

  script:
    """
    host=\$(hostname)
    echo \${host}

    bedtools \
    groupby \
    -i ${bed} \
    ${groupby_params} \
    ${params.bedtoolsgroupby.optional} \
    > ${sample_id}.${feature_id}.groupby.bed
    """
}

process groupbyAll {
  tag {"bedtoolsGroupbyAll ${sample_id}"}
  label 'bedtoolsGroupbyAll'
  label 'bedtools_2_30_0_groupbyAll'
  shell = ['/bin/bash', '-euo', 'pipefail']
  container = 'quay.io/biocontainers/bedtools:2.30.0--h468198e_3'

  input:
    tuple( val(donor_id), val(sample_id), path(bed))

  output:
    tuple( val(donor_id), val(sample_id), path("${sample_id}.features.groupby.bed"), emit: feature_groupby_bed)

  script:
    """
    host=\$(hostname)
    echo \${host}

    bedtools \
    groupby \
    -i ${bed} \
    ${params.bedtoolsgroupbyall.optional} \
    > ${sample_id}.features.groupby.bed
    """
}
