process createSnvRfTable {
  tag {"createSnvRfTable ${sample_id}"}
  label 'createSnvRfTable'
  shell = ['/bin/bash', '-euo', 'pipefail']
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://vanboxtelbioinformatics/ptato_r:1.3.3':
        'europe-west4-docker.pkg.dev/pmc-gcp-box-d-pip-development/pipeline-containers/ptato_r@sha256:c2396d1d9c217123444f4fed846fba38bb5dade14833b244a33bb9806577b059' }"

  input:
    tuple( val(donor_id), val(sample_id), path(ab_table), path(bed) )

  output:
    tuple( val(donor_id), val(sample_id), path("${sample_id}.rftable.rds"), emit: rf_table)

  script:
    """
    host=\$(hostname)
    echo \${host}

    R --slave --file=/scripts/R/create_snv_rf_table.R --args ${ab_table} ${bed} ${donor_id} ${sample_id}.rftable.rds
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

    R --slave --file=/scripts/R/create_indel_rf_table.R --args ${ab_table} ${bed} ${donor_id} ${sample_id}.rftable.rds
    """
}
