process closest {
  tag {"bedtoolsClosest ${sample_id} ${feature_id}"}
  label 'bedtoolsClosest'
  label 'bedtools_2_30_0_closest'
  shell = ['/bin/bash', '-euo', 'pipefail']
  container = 'quay.io/biocontainers/bedtools:2.30.0--h468198e_3'

  input:
    tuple( val(donor_id), val(sample_id), path(bed), val(feature_id), path(feature_bed), val(merge_params))

  output:
    tuple( val(donor_id), val(sample_id), val(feature_id), path("${sample_id}.${feature_id}.bed"), val(merge_params), emit: feature_bed)

  script:
    """
    host=\$(hostname)
    echo \${host}

    bedtools \
    closest \
    -a ${bed} \
    -b ${feature_bed} \
    ${params.bedtoolsclosest.optional} \
    > ${sample_id}.${feature_id}.bed
    """
}
