process createABtable {
  tag {"createABtable ${sample_id} ${chrom}"}
  label 'createABtable'
  shell = ['/bin/bash', '-euo', 'pipefail']

  input:
    tuple( val(donor_id), val(sample_id), val(chrom), path(phased_vcf), path(phased_tbi), path(germline_vcf), path(germline_tbi), path(somatic_vcf), path(somatic_tbi), val(bulk_name) )

  output:
    tuple( val(donor_id), val(sample_id), path("${sample_id}_${chrom}.abtable.txt"), emit: ab_table)

  script:
    """
    host=\$(hostname)
    echo \${host}

    R --slave --file=${baseDir}/scripts/R/ABscript.R --args ${somatic_vcf} ${germline_vcf} ${phased_vcf} ${chrom} ${bulk_name} ${sample_id}_${chrom}.abtable.txt ${params.ref_genome}
    """
}

process mergeABtable {
  tag {"mergeABtable ${sample_id}"}
  label 'mergeABtable'
  shell = ['/bin/bash', '-euo', 'pipefail']

  input:
    tuple( val(donor_id), val(sample_id), path(ab_tables) )

  output:
    tuple( val(donor_id), val(sample_id), path("${sample_id}.abtable.txt"), emit: ab_table )

  script:
    """
    header=true
    for ab_table in ${ab_tables}; do
      if \$header; then
        cat \$ab_table > ${sample_id}.abtable.txt
      else
        cat \$ab_table | tail -n+2 >> ${sample_id}.abtable.txt
      fi
      header=false
    done
    """
}
