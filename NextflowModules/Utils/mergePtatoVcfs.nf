process mergePtatoVcfs {
  tag {"mergePtatoVcfs ${donor_id}"}
  label 'mergePtatoVcfs'
  shell = ['/bin/bash', '-euo', 'pipefail']

  input:
    tuple( val(donor_id), val(sample_id), path(input_vcf), path(input_tbi), val(ptato_snvs_sample_ids), path(ptato_snvs_vcfs), path(ptato_snvs_tbis), val(ptato_indels_sample_ids), path(ptato_indels_vcfs), path(ptato_indels_tbis) )

  output:
    tuple( val(donor_id), val(donor_id), path("${donor_id}.ptato.merged.vcf"), emit: ptato_intersect_vcf)

  script:
    ptato_vcfs = ptato_snvs_vcfs.join(',')+","+ptato_indels_vcfs.join(',')
    """
    host=\$(hostname)
    echo \${host}

    R --slave --file=${baseDir}/scripts/R/merge_ptato_vcfs.R --args ${input_vcf} ${ptato_vcfs} ${donor_id}.ptato.merged.vcf
    """
}
