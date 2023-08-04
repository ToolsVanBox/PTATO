process QC_PtatoPlus {
  tag {"QC_PtatoPlus"}
  label 'QC_PtatoPlus'
  shell = ['/bin/bash', '-euo', 'pipefail']

  input:
    tuple( val(donor_id), val(sample_ids), path(PTATO_files), path(PTATO_filt_files) )

  output:
    tuple( val(donor_id), path("${donor_id}.postPTATO_QC.pdf"), path("${donor_id}.postPTATO_QC.txt") )

  script:
    input_args_1 = sample_ids.join(',')
    input_args_2 = PTATO_files.join(',')
    input_args_3 = PTATO_filt_files.join(',')
    """
    host=\$(hostname)
    echo \${host}

    R --slave --file=${baseDir}/scripts/R/PTATO_Report.R --args ${input_args_1} ${input_args_2} ${input_args_3} ${donor_id}.postPTATO_QC
    """

}