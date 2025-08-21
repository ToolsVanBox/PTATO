process ptatoFilter {
  tag {"ptatoFilter ${sample_id}"}
  label 'ptatoFilter'
  shell = ['/bin/bash', '-euo', 'pipefail']
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?  
        'docker://vanboxtelbioinformatics/ptato_r:1.3.3':
        'europe-west4-docker.pkg.dev/pmc-gcp-box-d-pip-development/pipeline-containers/ptato_rsha256:01608f437f46324b7ef4c6ae3ce9fe06da8d160ab6d46a44b04342e979f13be6' }"

  input:
    tuple( val(donor_id), val(sample_id), path(ptaprob_cutoff), path(ptato_vcf), path(ptato_tbi), path(walker_vcf), path(walker_tbi))

  output:
    tuple( val(donor_id), val(sample_id), path("${sample_id}.snvs.ptato.filtered.vcf"), emit: ptatofilter_out)

  script:
    """
    host=\$(hostname)
    echo \${host}

    R --slave --file=/scripts/R/ptatoFilter.R --args ${ptato_vcf} ${walker_vcf} ${ptaprob_cutoff} ${sample_id}.snvs.ptato.filtered.vcf
    """
}

process ptatoIndelFilter {
  tag {"ptatoIndelFilter ${sample_id}"}
  label 'ptatoIndelFilter'
  shell = ['/bin/bash', '-euo', 'pipefail']

  input:
    tuple( val(donor_id), val(sample_id), path(vcf), path(tbi) )

  output:
    tuple( val(donor_id), val(sample_id), path("${sample_id}.indels.ptato.filtered.vcf"), emit: ptatofilter_out)

  script:
    """
    host=\$(hostname)
    echo \${host}

    R --slave --file=/scripts/R/ptatoIndelFilter.R --args ${vcf} ${sample_id}.indels.ptato.filtered.vcf
    """
}
