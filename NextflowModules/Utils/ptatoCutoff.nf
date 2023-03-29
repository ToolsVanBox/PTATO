process ptatoCutoff {
  tag {"ptatoCutoff ${sample_id}"}
  label 'ptatoCutoff'
  shell = ['/bin/bash', '-euo', 'pipefail']

  input:
    tuple( val(donor_id), val(sample_id), path(ptato_vcf), path(ptato_tbi), path(walker_vcf), path(walker_tbi))

  output:
    tuple( val(donor_id), val(sample_id), path("${sample_id}.ptatotable.txt"), path("${sample_id}.ptaprobcutoff.txt"), emit: ptatofilter_out)

  script:
    """
    host=\$(hostname)
    echo \${host}

    R --slave --file=${baseDir}/scripts/R/ptatoCutoff.R --args ${ptato_vcf} ${walker_vcf} ${params.ref_genome} ${params.ptatocutoff.optional} ${sample_id}.ptatotable.txt > ${sample_id}.ptaprobcutoff.txt
    """
}
