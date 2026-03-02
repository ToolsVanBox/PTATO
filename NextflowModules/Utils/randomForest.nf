process train_snv_rf {
  tag {"train_snv_rf"}
  label 'train_snv_rf'
  shell = ['/bin/bash', '-euo', 'pipefail']
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?  
        'docker://vanboxtelbioinformatics/ptato_r:1.3.3':
        params.artifact_registry_path + '/ptato_r:1.3.3' }"
//        'europe-west4-docker.pkg.dev/pmc-gcp-box-d-pip-development/pipeline-containers/ptato_r@sha256:c2396d1d9c217123444f4fed846fba38bb5dade14833b244a33bb9806577b059' }"

  input:
    tuple( val(label_1), path(rf_table_1), val(label_2), path(rf_table_2) )

  output:
    tuple( path("randomforest*confusion.txt"), path("randomforest*importance.txt"), path("randomforest*.Rdata"), path("randomforest*.rds"), emit: random_forest_file )

  script:
    input_args_1 = rf_table_1.join(',')
    input_args_2 = rf_table_2.join(',')
    """
    host=\$(hostname)
    echo \${host}

    R --slave --file=/scripts/R/train_randomforest.R --args ${label_1} ${input_args_1} ${label_2} ${input_args_2} ${params.train.version}
    """
}

process train_indel_rf {
  tag {"train_indel_rf"}
  label 'train_indel_rf'
  shell = ['/bin/bash', '-euo', 'pipefail']
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?  
        'docker://vanboxtelbioinformatics/ptato_r:1.3.3':
        params.artifact_registry_path + '/ptato_r:1.3.3' }"
//        'europe-west4-docker.pkg.dev/pmc-gcp-box-d-pip-development/pipeline-containers/ptato_r@sha256:c2396d1d9c217123444f4fed846fba38bb5dade14833b244a33bb9806577b059' }"

  input:
    tuple( val(label_1), path(rf_table_1), val(label_2), path(rf_table_2) )

  output:
    tuple( path("randomforest*confusion.txt"), path("randomforest*importance.txt"), path("randomforest*.Rdata"), path("randomforest*.rds"), emit: random_forest_file )

  script:
    input_args_1 = rf_table_1.join(',')
    input_args_2 = rf_table_2.join(',')
    """
    host=\$(hostname)
    echo \${host}

    R --slave --file=/scripts/R/train_randomforest.R --args ${label_1} ${input_args_1} ${label_2} ${input_args_2} ${params.train.version}
    """
}

process test_snv_rf {
  tag {"test_snv_rf ${sample_id}"}
  label 'test_snv_rf'
  shell = ['/bin/bash', '-euo', 'pipefail']
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?  
        'docker://vanboxtelbioinformatics/ptato_r:1.3.3':
        params.artifact_registry_path + '/ptato_r:1.3.3' }"
//        'europe-west4-docker.pkg.dev/pmc-gcp-box-d-pip-development/pipeline-containers/ptato_r@sha256:c2396d1d9c217123444f4fed846fba38bb5dade14833b244a33bb9806577b059' }"

  input:
    tuple( val(donor_id), val(sample_id), path(somatic_vcf), path(somatic_tbi), path(rf_table) )
    path( rf_rds )
  output:
    tuple( val(donor_id), val(sample_id), path("${sample_id}.snvs.ptato.vcf"), emit: ptato_vcf )

  script:
    """
    host=\$(hostname)
    echo \${host}

    R --slave --file=/scripts/R/test_randomforest.R --args ${rf_rds} ${rf_table} ${somatic_vcf} ${sample_id}.snvs.ptato.vcf
    """
}

process test_indel_rf {
  tag {"test_indel_rf ${sample_id}"}
  label 'test_indel_rf'
  shell = ['/bin/bash', '-euo', 'pipefail']
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?  
        'docker://vanboxtelbioinformatics/ptato_r:1.3.3':
        params.artifact_registry_path + '/ptato_r:1.3.3' }"
//        'europe-west4-docker.pkg.dev/pmc-gcp-box-d-pip-development/pipeline-containers/ptato_r@sha256:c2396d1d9c217123444f4fed846fba38bb5dade14833b244a33bb9806577b059' }"

  input:
    tuple( val(donor_id), val(sample_id), path(somatic_vcf), path(somatic_tbi), path(rf_table) )
    path( rf_rds )
  output:
    tuple( val(donor_id), val(sample_id), path("${sample_id}.indels.ptato.vcf"), emit: ptato_vcf )

  script:
    """
    host=\$(hostname)
    echo \${host}

    R --slave --file=/scripts/R/test_randomforest.R --args ${rf_rds} ${rf_table} ${somatic_vcf} ${sample_id}.indels.ptato.vcf
    """
}
