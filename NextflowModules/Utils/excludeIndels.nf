process excludeIndels {
  tag {"excludeIndels ${sample_id}"}
  label 'excludeIndels'
  shell = ['/bin/bash', '-euo', 'pipefail']
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://vanboxtelbioinformatics/ptato_r:1.3.3':
        'europe-west4-docker.pkg.dev/pmc-gcp-box-d-pip-development/pipeline-containers/ptato_rsha256:01608f437f46324b7ef4c6ae3ce9fe06da8d160ab6d46a44b04342e979f13be6' }"

  input:
    tuple( val(donor_id), val(sample_id), path(indels_vcf), path(indels_tbi), path(context_bed), path(exclude_vcf) )

  output:
  tuple( val(donor_id), val(sample_id), path("${sample_id}.indels.ptato.vcf"), emit: ptato_vcf )


  script:
    """
    host=\$(hostname)
    echo \${host}

    R --slave --file=/scripts/R/excludeIndels.R --args ${indels_vcf} ${context_bed} ${exclude_vcf} ${sample_id}.indels.ptato.vcf
    """
}
