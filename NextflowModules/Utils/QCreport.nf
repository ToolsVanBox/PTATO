process QCreport {
  tag {"QCreport"}
  label 'QCreport'
  shell = ['/bin/bash', '-euo', 'pipefail']
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?  
        'docker://vanboxtelbioinformatics/ptato_r:1.3.3':
        'europe-west4-docker.pkg.dev/pmc-gcp-box-d-pip-development/pipeline-containers/ptato_r@sha256:c2396d1d9c217123444f4fed846fba38bb5dade14833b244a33bb9806577b059' }"

  input:
    tuple( val(donor_id), val(sample_ids), path(insert_size_metrics_files), path(wgs_metrics_files) )

  output:
    tuple( val(donor_id), path("${donor_id}.qcreport.pdf"), path("${donor_id}.qcreport.txt"), emit: qc_report_pdf )

  script:
    input_args_1 = sample_ids.join(',')
    input_args_2 = insert_size_metrics_files.join(',')
    input_args_3 = wgs_metrics_files.join(',')
    """
    host=\$(hostname)
    echo \${host}

    R --slave --file=/scripts/R/PTA_QC_report.R --args ${input_args_1} ${input_args_2} ${input_args_3} ${donor_id}.qcreport.pdf
    """
}
