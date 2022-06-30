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

process intersectPON {
  tag {"bedtoolsIntersectPON ${sample_id}"}
  label 'bedtoolsIntersectPON'
  label 'bedtools_2_30_0_intersectPON'
  shell = ['/bin/bash', '-euo', 'pipefail']
  container = 'quay.io/biocontainers/bedtools:2.30.0--h468198e_3'

  input:
    tuple( val(donor_id), val(sample_id), path(vcf), path(tbi), path(pon) )

  output:
    tuple( val(donor_id), val(sample_id), path("${sample_id}.pon.filtered.vcf"), emit: pon_filtered_vcf)

  script:
    """
    host=\$(hostname)
    echo \${host}

    bedtools \
    intersect \
    -a ${vcf} \
    -b ${pon} \
    ${params.bedtoolsintersectpon.optional} \
    > ${sample_id}.pon.filtered.vcf
    """
}

process intersectPTATO {
  tag {"bedtoolsIntersectPTATO ${donor_id}"}
  label 'bedtoolsIntersectPTATO'
  label 'bedtools_2_30_0_intersectPTATO'
  shell = ['/bin/bash', '-euo', 'pipefail']
  container = 'quay.io/biocontainers/bedtools:2.30.0--h468198e_3'

  input:
    tuple( val(donor_id), val(sample_id), path(input_vcf), path(input_tbi), val(ptato_snvs_sample_ids), path(ptato_snvs_vcfs), path(ptato_snvs_tbis), val(ptato_indels_sample_ids), path(ptato_indels_vcfs), path(ptato_indels_tbis) )

  output:
    tuple( val(donor_id), val(donor_id), path("${donor_id}.ptato.intersect.vcf"), emit: ptato_intersect_vcf)

  script:
    b1 = ptato_snvs_vcfs ? ' -b ' + ptato_snvs_vcfs.join(' -b ') : ''
    b2 = ptato_indels_vcfs ? ' -b ' + ptato_indels_vcfs.join(' -b ') : ''
    """
    host=\$(hostname)
    echo \${host}

    bedtools \
    intersect \
    -a ${input_vcf} \
    ${b1}\
    ${b2}\
    ${params.bedtoolsintersectptato.optional} \
    > ${donor_id}.ptato.intersect.vcf
    """
}
