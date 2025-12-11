process AutosomalCallableLoci {
  tag {"AutosomalCallableLoci"}
  label 'AutosomalCallableLoci'
  shell = ['/bin/bash', '-euo', 'pipefail']
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://vanboxtelbioinformatics/ptato_r:1.3.3':
        'europe-west4-docker.pkg.dev/pmc-gcp-box-d-pip-development/pipeline-containers/ptato_r@sha256:c2396d1d9c217123444f4fed846fba38bb5dade14833b244a33bb9806577b059' }"


  input:
    tuple( val(donor_id), val(sample_id), path(callable_bed_file), path(callable_txt_file) )
  output:
    tuple( val(donor_id), val(sample_id), path("${sample_id}.callableloci.autosomal.txt"), emit: autosomal_callable_file )
  script:
    """
    host=\$(hostname)
    echo \${host}

    R --slave --file=/scripts/R/AutosomalCallableLoci.R --args ${callable_bed_file} ${sample_id}.callableloci.autosomal.txt
    """
}