process ptatoFilter {
  tag {"ptatoFilter ${sample_id}"}
  label 'ptatoFilter'
  shell = ['/bin/bash', '-euo', 'pipefail']

  input:
    tuple( val(donor_id), val(sample_id), path(ptato_vcf), path(ptato_tbi), path(walker_vcf), path(walker_tbi))

  output:
    tuple( val(donor_id), val(sample_id), path("${sample_id}.ptato.filtered.vcf"), path("${sample_id}.ptatotable.txt"), emit: ptatofilter_out)

  script:
    """
    host=\$(hostname)
    echo \${host}

    R --slave --file=${baseDir}/scripts/R/ptatoFilter.R --args ${ptato_vcf} ${walker_vcf} ${sample_id}.ptato.filtered.vcf ${sample_id}.ptatotable.txt
    """
}
