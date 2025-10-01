process getContext {
  tag {"getContext ${sample_id}"}
  label 'getContext'
  shell = ['/bin/bash', '-euo', 'pipefail']
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?  
        'docker://vanboxtelbioinformatics/ptato_r:1.3.3':
        'europe-west4-docker.pkg.dev/pmc-gcp-box-d-pip-development/pipeline-containers/ptato_r@sha256:c2396d1d9c217123444f4fed846fba38bb5dade14833b244a33bb9806577b059' }"


  input:
    tuple( val(donor_id), val(sample_id), path(vcf), path(tbi) )

  output:
    tuple( val(donor_id), val(sample_id), path("${sample_id}.context.bed"), emit: bed)

  script:
    """
    host=\$(hostname)
    echo \${host}

    R --slave --file=/scripts/R/get_context.R --args ${vcf} ${sample_id}.context.bed ${params.ref_genome}
    """
}
