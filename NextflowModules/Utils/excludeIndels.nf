process excludeIndels {
  tag {"excludeIndels ${sample_id}"}
  label 'excludeIndels'
  shell = ['/bin/bash', '-euo', 'pipefail']

  input:
    tuple( val(donor_id), val(sample_id), path(indels_vcf), path(indels_tbi), path(context_bed), path(exclude_vcf) )

  output:
  tuple( val(donor_id), val(sample_id), path("${sample_id}.indels.ptato.vcf"), emit: ptato_vcf )


  script:
    """
    host=\$(hostname)
    echo \${host}

    R --slave --file=${baseDir}/scripts/R/excludeIndels.R --args ${indels_vcf} ${context_bed} ${exclude_vcf} ${sample_id}.indels.ptato.vcf
    """
}
