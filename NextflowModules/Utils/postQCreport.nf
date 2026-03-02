process postQCreport {
  tag {"postQCreport"}
  label 'postQCreport'
  shell = ['/bin/bash', '-euo', 'pipefail']
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?  
        'docker://vanboxtelbioinformatics/ptato_r:1.3.3':
        params.artifact_registry_path + '/ptato_r:1.3.3' }"
//        'europe-west4-docker.pkg.dev/pmc-gcp-box-d-pip-development/pipeline-containers/ptato_r@sha256:c2396d1d9c217123444f4fed846fba38bb5dade14833b244a33bb9806577b059' }"

  input:
    tuple( val(donor_id), val(sample_ids), path(ptato_vcf), path(ptato_tbi), path(ptato_filt_vcf), path(ptato_filt_tbi), path(ptato_table), path(walker_vcf), path(walker_tbi), path(autosomal_callable_files) )

  output:
    tuple( val(donor_id), val(sample_ids), path("${donor_id}.postqcreport.pdf"), path("${donor_id}.postqcreport.txt"), emit: postqc_report_pdf )

  script:
    input_args_1 = sample_ids.join(',')
    input_args_2 = ptato_vcf.join(',')
    input_args_3 = ptato_filt_vcf.join(',')
    input_args_4 = walker_vcf.join(',')
    input_args_5 = ptato_table.join(',')
    input_args_6 = autosomal_callable_files.join(',')

    """
    host=\$(hostname)
    echo \${host}

    R --slave --file=/scripts/R/PostPTATO_QC_report.R --args ${input_args_1} ${input_args_2} ${input_args_3} ${input_args_4} ${input_args_5} ${input_args_6} ${donor_id}.postqcreport
    """
}