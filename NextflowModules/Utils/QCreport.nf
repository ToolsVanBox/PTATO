process QCreport {
  tag {"QCreport"}
  label 'QCreport'
  shell = ['/bin/bash', '-euo', 'pipefail']

  input:
    tuple( val(donor_id), val(sample_ids), path(insert_size_metrics_files), path(wgs_metrics_files) )

  output:
    tuple( val(donor_id), path("${donor_id}.qcreport.pdf"), emit: qc_report_pdf )

  script:
    input_args_1 = sample_ids.join(',')
    input_args_2 = insert_size_metrics_files.join(',')
    input_args_3 = wgs_metrics_files.join(',')
    """
    host=\$(hostname)
    echo \${host}

    R --slave --file=${baseDir}/scripts/R/PTA_QC_report.R --args ${input_args_1} ${input_args_2} ${input_args_3} ${donor_id}.qcreport.pdf
    """
}
