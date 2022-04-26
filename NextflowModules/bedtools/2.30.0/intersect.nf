process intersect {
  tag {"bedtoolsIntersect ${sample_id} ${feature_id}"}
  label 'bedtoolsIntersect'
  label 'bedtools_2_30_0_intersect'
  shell = ['/bin/bash', '-euo', 'pipefail']
  container = 'quay.io/biocontainers/bedtools:2.30.0--h468198e_3'

  input:
    tuple( val(donor_id), val(sample_id), path(bed), val(feature_id), path(feature_bed), val(merge_params) )

  output:
    tuple( val(donor_id), val(sample_id), val(feature_id), path("${sample_id}.${feature_id}.bed"), val(merge_params), emit: feature_bed)

  script:
    """
    host=\$(hostname)
    echo \${host}

    bedtools \
    intersect \
    -a ${bed} \
    -b ${feature_bed} \
    ${params.bedtoolsintersect.optional} \
    > ${sample_id}.${feature_id}.bed
    """
}

process intersectAll {
  tag {"bedtoolsIntersectAll ${sample_id}"}
  label 'bedtoolsIntersectAll'
  label 'bedtools_2_30_0_intersectAll'
  shell = ['/bin/bash', '-euo', 'pipefail']
  container = 'quay.io/biocontainers/bedtools:2.30.0--h468198e_3'

  input:
    tuple( val(donor_id), val(sample_id), path(bed), val(closest_feature_beds), val(intersect_feature_beds) )

  output:
    tuple( val(donor_id), val(sample_id), path("${sample_id}.features.bed"), emit: features_bed)

  script:
    b1 = closest_feature_beds ? ' -b ' + closest_feature_beds.join(' -b ') : ''
    b2 = intersect_feature_beds ? ' -b ' + intersect_feature_beds.join(' -b ') : ''
    """
    host=\$(hostname)
    echo \${host}
    names=''
    regex=".+\\.([A-Z]+).merged.bed"
    for feature in ${b1} ${b2}; do
      if [[ \$feature =~ \$regex ]]; then
        name="\${BASH_REMATCH[1]}"
        names=\$names' '\$name
      fi
    done

    bedtools \
    intersect \
    -a ${bed} \
    ${b1} \
    ${b2} \
    -names \${names} \
    ${params.bedtoolsintersectall.optional} \
    > ${sample_id}.features.bed
    """
}
