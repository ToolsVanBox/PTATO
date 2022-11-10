process ptatoFilter {
  tag {"ptatoFilter ${sample_id}"}
  label 'ptatoFilter'
  shell = ['/bin/bash', '-euo', 'pipefail']

  input:
    tuple( val(donor_id), val(sample_id), path(ptaprob_cutoff), path(ptato_vcf), path(ptato_tbi), path(walker_vcf), path(walker_tbi))

  output:
    tuple( val(donor_id), val(sample_id), path("${sample_id}.snvs.ptato.filtered.vcf"), emit: ptatofilter_out)

  script:
    """
    host=\$(hostname)
    echo \${host}

    R --slave --file=${baseDir}/scripts/R/ptatoFilter.R --args ${ptato_vcf} ${walker_vcf} ${ptaprob_cutoff} ${sample_id}.snvs.ptato.filtered.vcf
    """
}

process ptatoIndelFilter {
  tag {"ptatoIndelFilter ${sample_id}"}
  label 'ptatoIndelFilter'
  shell = ['/bin/bash', '-euo', 'pipefail']

  input:
    tuple( val(donor_id), val(sample_id), path(vcf), path(tbi) )

  output:
    tuple( val(donor_id), val(sample_id), path("${sample_id}.indels.ptato.filtered.vcf"), emit: ptatofilter_out)

  script:
    """
    host=\$(hostname)
    echo \${host}

    R --slave --file=${baseDir}/scripts/R/ptatoIndelFilter.R --args ${vcf} ${sample_id}.indels.ptato.filtered.vcf
    """
}
