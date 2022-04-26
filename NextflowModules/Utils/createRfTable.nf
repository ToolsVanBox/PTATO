process createSnvRfTable {
  tag {"createSnvRfTable ${sample_id}"}
  label 'createSnvRfTable'
  shell = ['/bin/bash', '-euo', 'pipefail']

  input:
    tuple( val(donor_id), val(sample_id), path(ab_table), path(bed) )

  output:
    tuple( val(donor_id), val(sample_id), path("${sample_id}.rftable.rds"), emit: rf_table)

  script:
    """
    host=\$(hostname)
    echo \${host}

    cat ${baseDir}/scripts/R/create_snv_rf_table.R | R --slave --args ${ab_table} ${bed} ${donor_id} ${sample_id}.rftable.rds
    """
}

process createIndelRfTable {
  tag {"createSnvRfTable ${sample_id}"}
  label 'createSnvRfTable'
  shell = ['/bin/bash', '-euo', 'pipefail']

  input:
    tuple( val(donor_id), val(sample_id), path(ab_table), path(bed) )

  output:
    tuple( val(donor_id), val(sample_id), path("${sample_id}.rftable.rds"), emit: rf_table)

  script:
    """
    host=\$(hostname)
    echo \${host}

    cat ${baseDir}/scripts/R/create_indel_rf_table.R | R --slave --args ${ab_table} ${bed} ${donor_id} ${sample_id}.rftable.rds
    """
}
