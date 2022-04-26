process train_snv_rf {
  tag {"train_snv_rf"}
  label 'train_snv_rf'
  shell = ['/bin/bash', '-euo', 'pipefail']

  input:
    tuple( val(label_1), path(rf_table_1), val(label_2), path(rf_table_2) )

  output:
    tuple( path("randomforest_v1.0.0_confusion.txt"), path("randomforest_v1.0.0_importance.txt"), path("randomforest_v1.0.0.Rdata"), path("randomforest_v1.0.0.rds"), emit: random_forest_file )

  script:
    input_args_1 = rf_table_1.join(',')
    input_args_2 = rf_table_2.join(',')
    """
    host=\$(hostname)
    echo \${host}

    cat ${baseDir}/scripts/R/train_randomforest.R | R --slave --args ${label_1} ${input_args_1} ${label_2} ${input_args_2}
    """
}

process train_indel_rf {
  tag {"train_indel_rf"}
  label 'train_indel_rf'
  shell = ['/bin/bash', '-euo', 'pipefail']

  input:
    tuple( val(label_1), path(rf_table_1), val(label_2), path(rf_table_2) )

  output:
    tuple( path("randomforest_v1.0.0_confusion.txt"), path("randomforest_v1.0.0_importance.txt"), path("randomforest_v1.0.0.Rdata"), path("randomforest_v1.0.0.rds"), emit: random_forest_file )

  script:
    input_args_1 = rf_table_1.join(',')
    input_args_2 = rf_table_2.join(',')
    """
    host=\$(hostname)
    echo \${host}

    cat ${baseDir}/scripts/R/train_randomforest.R | R --slave --args ${label_1} ${input_args_1} ${label_2} ${input_args_2}
    """
}

process test_snv_rf {
  tag {"test_snv_rf ${sample_id}"}
  label 'test_snv_rf'
  shell = ['/bin/bash', '-euo', 'pipefail']

  input:
    tuple( val(donor_id), val(sample_id), path(somatic_vcf), path(somatic_tbi), path(rf_table) )

  output:
    tuple( val(donor_id), val(sample_id), path("${sample_id}.ptato.vcf"), emit: pap_vcf )

  script:
    """
    host=\$(hostname)
    echo \${host}

    cat ${baseDir}/scripts/R/test_randomforest.R | R --slave --args ${params.snvs.rf_rds} ${rf_table} ${somatic_vcf} ${sample_id}.ptato.vcf
    """
}

process test_indel_rf {
  tag {"test_indel_rf ${sample_id}"}
  label 'test_indel_rf'
  shell = ['/bin/bash', '-euo', 'pipefail']

  input:
    tuple( val(donor_id), val(sample_id), path(somatic_vcf), path(somatic_tbi), path(rf_table) )

  output:
    tuple( val(donor_id), val(sample_id), path("${sample_id}.ptato.vcf"), emit: pap_vcf )

  script:
    """
    host=\$(hostname)
    echo \${host}

    cat ${baseDir}/scripts/R/test_randomforest.R | R --slave --args ${params.indels.rf_rds} ${rf_table} ${somatic_vcf} ${sample_id}.ptato.vcf
    """
}
